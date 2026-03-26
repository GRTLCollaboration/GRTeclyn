/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SUPPORTEDWORMHOLEINITIALDATA_HPP_
#define SUPPORTEDWORMHOLEINITIALDATA_HPP_

#include "CCZ4StateVariables.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "Tensor.hpp"
#include "VarsTools.hpp"
#include "simd.hpp"

#include <cmath>

//! Initial data for a self-consistent single-throat Ellis-Bronnikov wormhole.
//!
//! Uses the conformally-flat isotropic form in Cartesian coordinates:
//!
//!   psi = sqrt(1 + b0^2 / (4 r^2)),   chi = psi^{-4},   h_ij = delta_ij.
//!
//! Equivalently: chi = (1 + b0^2/(4r^2))^{-2} = (4r^2/(4r^2 + b0^2))^2.
//!
//! The evolution system is the fully coupled Einstein-Klein-Gordon system
//! with a phantom (ghost) scalar field. This initial slice (with Pi=0,
//! K_ij=0) satisfies the momentum constraint exactly. The Hamiltonian
//! constraint has only a small spin-0 residual from the scalar perturbation.
class SupportedWormholeInitialData
{
  public:
    struct params_t
    {
        // Initial lapse selector:
        // 0 = alpha = 1
        // 1 = alpha = sqrt(chi) (pre-collapsed)
        // 2 = alpha = 1 - 3*log(chi) (N=2, 1+log-inspired using gamma = chi^{-3})
        int initial_lapse_type;

        // Grid center used for index->physical coordinate mapping
        std::array<double, AMREX_SPACEDIM> grid_center;

        // Mouth parameter
        double b0; // Throat radius
        std::array<double, AMREX_SPACEDIM> centerA;

        // Scalar field profile perturbation at t=0:
        //   phi -> phi_EB(r) + (A_0 + A_phi * Y20(theta)) * exp(-r^2/sigma^2)
        // With Pi=0 and K_ij=0, the momentum constraint is satisfied exactly.
        // Only the Hamiltonian constraint picks up a small spin-0 residual
        // from the changed gradient energy, which does not contaminate Psi4.
        double phi_monopole_amplitude;
        double phi_perturbation_amplitude;
        double phi_perturbation_width;

        double phantom_mass;
        double support_strength;
    };

    SupportedWormholeInitialData(params_t a_params, double a_dx)
        : m_params(a_params), m_dx(a_dx)
    {
    }

    template <class data_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    compute(int i, int j, int k, amrex::Array4<data_t> cell) const
    {
        amrex::IntVect grid_index(i, j, k);
        Coordinates coords(grid_index, m_dx, m_params.grid_center);

        const data_t x = coords.x;
        const data_t y = coords.y;
        const data_t z = coords.z;

        const double b0 = m_params.b0;
        const double b0_sq = b0 * b0;

        // Distances to throat
        const data_t dxA = x - (data_t)m_params.centerA[0];
        const data_t dyA = y - (data_t)m_params.centerA[1];
        const data_t dzA = z - (data_t)m_params.centerA[2];
        const data_t rA2 = dxA * dxA + dyA * dyA + dzA * dzA;

        // Regularisation to avoid division by zero at origin
        const data_t eps2 = (data_t)1.0e-24;
        const data_t rA2_reg = simd_max(rA2, eps2);

        data_t chi = 1.0;
        data_t h11 = 1.0, h12 = 0.0, h13 = 0.0, h22 = 1.0, h23 = 0.0, h33 = 1.0;

        // Isotropic conformally-flat Ellis–Bronnikov form:
        //   psi^2 = 1 + b0^2/(4r^2),  chi = psi^{-4} = (1 + b0^2/(4r^2))^{-2}
        const data_t termA = (data_t)b0_sq / (4.0 * rA2_reg);
        const data_t psi_sq = 1.0 + termA;
        chi = 1.0 / (psi_sq * psi_sq);

        // Floors (avoid NaNs in evolution)
        if (chi < (data_t)1.0e-10) chi = (data_t)1.0e-10;

        // Initialize scalar field phi and Pi
        data_t phi = 0.0;
        data_t Pi = 0.0;
        
        // phi(r) = (1/sqrt(4*pi)) * arctan( (r - b0^2/(4r)) / b0 )
        const data_t rA = sqrt(rA2_reg);
        const data_t argA = (rA - (data_t)b0_sq / (4.0 * rA)) / (data_t)b0;
        phi = (data_t)(1.0 / sqrt(4.0 * M_PI)) * atan(argA);

        // Scalar field profile perturbation at t=0:
        //   phi += (A_0 + A_phi * Y20(theta)) * exp(-r^2 / sigma^2)
        // A_0 = monopolar amplitude (drives fast symmetric collapse)
        // A_phi = quadrupolar amplitude (breaks spherical symmetry -> GW signal)
        // Pi remains zero, so the momentum constraint is exactly satisfied.
        const double phi_mono = m_params.phi_monopole_amplitude;
        const double phi_amp = m_params.phi_perturbation_amplitude;
        const double phi_width = m_params.phi_perturbation_width;
        if ((phi_mono != 0.0 || phi_amp != 0.0) && phi_width > 0.0)
        {
            const data_t sig2 = (data_t)(phi_width * phi_width);
            const data_t eps2_ang = (data_t)1.0e-24;

            const data_t rA2_safe = simd_max(rA2, eps2_ang);
            const data_t cos_theta_A_sq = dzA * dzA / rA2_safe;
            const data_t Y20_A = (data_t)sqrt(5.0 / (16.0 * M_PI))
                               * (3.0 * cos_theta_A_sq - 1.0);
            const data_t dphi_A = ((data_t)phi_mono + (data_t)phi_amp * Y20_A)
                                  * exp(-rA2 / sig2);

            phi += dphi_A;
        }

        data_t lapse = 1.0;
        if (m_params.initial_lapse_type == 1)
        {
            lapse = sqrt(chi);
        }
        else if (m_params.initial_lapse_type == 2)
        {
            // For conformal variables with det(h)=1, det(gamma)=chi^{-3}.
            // Using alpha = 1 + ln(gamma^{N/2}) with N=2 gives:
            // alpha = 1 + ln(gamma) = 1 - 3 ln(chi).
            lapse = 1.0 - (data_t)3.0 * log(chi);
        }
        if (lapse < (data_t)1.0e-10) lapse = (data_t)1.0e-10;

        cell(i, j, k, c_chi) = chi;
        cell(i, j, k, c_h11) = h11;
        cell(i, j, k, c_h12) = h12;
        cell(i, j, k, c_h13) = h13;
        cell(i, j, k, c_h22) = h22;
        cell(i, j, k, c_h23) = h23;
        cell(i, j, k, c_h33) = h33;

        cell(i, j, k, c_K) = 0.0;
        cell(i, j, k, c_A11) = 0.0;
        cell(i, j, k, c_A12) = 0.0;
        cell(i, j, k, c_A13) = 0.0;
        cell(i, j, k, c_A22) = 0.0;
        cell(i, j, k, c_A23) = 0.0;
        cell(i, j, k, c_A33) = 0.0;

        cell(i, j, k, c_lapse)  = lapse;
        cell(i, j, k, c_shift1) = 0.0;
        cell(i, j, k, c_shift2) = 0.0;
        cell(i, j, k, c_shift3) = 0.0;
        cell(i, j, k, c_B1)     = 0.0;
        cell(i, j, k, c_B2)     = 0.0;
        cell(i, j, k, c_B3)     = 0.0;

        cell(i, j, k, c_Theta) = 0.0;
        
        cell(i, j, k, c_phi) = phi;
        cell(i, j, k, c_Pi) = Pi;
    }

  protected:
    params_t m_params;
    double m_dx;
};

#endif /* SUPPORTEDWORMHOLEINITIALDATA_HPP_ */