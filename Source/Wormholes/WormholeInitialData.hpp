/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef WORMHOLEINITIALDATA_HPP_
#define WORMHOLEINITIALDATA_HPP_

#include "CCZ4StateVariables.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "Tensor.hpp"
#include "VarsTools.hpp"
#include "simd.hpp"

#include <cmath>

//! Initial data for an "unsupported wormhole" experiment.
//!
//! This class provides a conformally-flat isotropic *two-mouth* ansatz via
//! superposition of Ellis–Bronnikov potentials:
//!
//!   psi(x) = 1 + bA^2/(4|x-cA|^2) + bB^2/(4|x-cB|^2),
//!   chi = psi^{-4},  h_ij = delta_ij.
//!
//! The evolution system is vacuum CCZ4. Therefore, this initial slice is not in
//! general a vacuum-constraint solution and will emit an early transient as it
//! relaxes.
class WormholeInitialData
{
  public:
    struct params_t
    {
        // Metric realisation selector:
        // 0 = two-mouth isotropic superposition (default, centers A/B)
        // 1 = single-throat Morris–Thorne / Ellis–Bronnikov (uses centerA only)
        int metric_type;

        // Initial lapse selector:
        // 0 = alpha = 1
        // 1 = alpha = sqrt(chi) (pre-collapsed)
        // 2 = alpha = 1 - 3*log(chi) (N=2, 1+log-inspired using gamma = chi^{-3})
        int initial_lapse_type;

        // Grid center used for index->physical coordinate mapping
        std::array<double, AMREX_SPACEDIM> grid_center;

        // Mouth parameters
        double throat_radius_A; // b_A
        double throat_radius_B; // b_B
        std::array<double, AMREX_SPACEDIM> centerA;
        std::array<double, AMREX_SPACEDIM> centerB;

        // Optional inward "kick" via a small negative K perturbation
        // applied as a sum of Gaussians centered at each mouth:
        //   K = k_amplitude * (exp(-rA^2/k_width^2) + exp(-rB^2/k_width^2))
        // Set k_amplitude = 0 to disable.
        double k_amplitude;
        double k_width;

        // Legacy/debug option: initialise a nontrivial Cartesian gamma_ij
        // (single-throat proper-distance metric centred at the origin).
        bool use_cartesian_gamma;
    };

    WormholeInitialData(params_t a_params, double a_dx)
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

        const double bA = m_params.throat_radius_A;
        const double bB = m_params.throat_radius_B;
        const double bA_sq = bA * bA;
        const double bB_sq = bB * bB;

        // Distances to mouths
        const data_t dxA = x - (data_t)m_params.centerA[0];
        const data_t dyA = y - (data_t)m_params.centerA[1];
        const data_t dzA = z - (data_t)m_params.centerA[2];
        const data_t rA2 = dxA * dxA + dyA * dyA + dzA * dzA;

        const data_t dxB = x - (data_t)m_params.centerB[0];
        const data_t dyB = y - (data_t)m_params.centerB[1];
        const data_t dzB = z - (data_t)m_params.centerB[2];
        const data_t rB2 = dxB * dxB + dyB * dyB + dzB * dzB;

        // Regularisation to avoid division by zero at a mouth center
        const data_t eps2 = (data_t)1.0e-24;
        const data_t rA2_reg = simd_max(rA2, eps2);
        const data_t rB2_reg = simd_max(rB2, eps2);

        data_t chi = 1.0;
        data_t h11 = 1.0, h12 = 0.0, h13 = 0.0, h22 = 1.0, h23 = 0.0, h33 = 1.0;

        if (m_params.use_cartesian_gamma)
        {
            // Legacy/debug: proper-distance single-throat metric at origin
            const data_t r0_2 = x * x + y * y + z * z;
            const data_t ell2 = r0_2 + (data_t)m_dx * (data_t)m_dx;
            const data_t ell  = sqrt(ell2);

            const data_t fac = (data_t)bA_sq / ell2;
            const data_t A   = 1.0 + fac;
            const data_t nx = x / ell;
            const data_t ny = y / ell;
            const data_t nz = z / ell;

            const data_t g11 = A - fac * nx * nx;
            const data_t g22 = A - fac * ny * ny;
            const data_t g33 = A - fac * nz * nz;
            const data_t g12 = -fac * nx * ny;
            const data_t g13 = -fac * nx * nz;
            const data_t g23 = -fac * ny * nz;

            // det(gamma)=A^2 (radial eigenvalue ~1) -> chi = A^{-2/3}
            chi = pow(A, (data_t)(-2.0 / 3.0));
            h11 = chi * g11;
            h22 = chi * g22;
            h33 = chi * g33;
            h12 = chi * g12;
            h13 = chi * g13;
            h23 = chi * g23;
        }
        else
        {
            // Isotropic conformally-flat Ellis–Bronnikov form:
            // - metric_type=0: two-mouth superposition
            // - metric_type=1: single-throat Morris–Thorne / Ellis–Bronnikov
            const data_t termA = (data_t)bA_sq / (4.0 * rA2_reg);
            const data_t termB = (m_params.metric_type == 1)
                                     ? (data_t)0.0
                                     : (data_t)bB_sq / (4.0 * rB2_reg);
            const data_t psi   = 1.0 + termA + termB;
            const data_t psi2  = psi * psi;
            const data_t psi4  = psi2 * psi2;
            chi = 1.0 / psi4;
        }

        // Floors (avoid NaNs in evolution)
        if (chi < (data_t)1.0e-10) chi = (data_t)1.0e-10;

        // Optional "kick": a small inward K perturbation at each mouth.
        data_t K = 0.0;
        if ((m_params.k_amplitude != 0.0) && (m_params.k_width > 0.0))
        {
            const data_t sig2 = (data_t)(m_params.k_width * m_params.k_width);
            if (m_params.metric_type == 1)
            {
                // Single throat: seed collapse at centerA only
                K = (data_t)m_params.k_amplitude * exp(-rA2 / sig2);
            }
            else
            {
                // Two mouths: kick both
                K = (data_t)m_params.k_amplitude *
                    (exp(-rA2 / sig2) + exp(-rB2 / sig2));
            }
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

        cell(i, j, k, c_K) = K;
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
    }

  protected:
    params_t m_params;
    double m_dx;
};

#endif /* WORMHOLEINITIALDATA_HPP_ */