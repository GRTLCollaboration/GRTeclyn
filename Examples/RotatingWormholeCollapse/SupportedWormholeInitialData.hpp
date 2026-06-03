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

class SupportedWormholeInitialData
{
  public:
    struct params_t
    {
        int initial_lapse_type;

        std::array<double, AMREX_SPACEDIM> grid_center;

        double b0;
        std::array<double, AMREX_SPACEDIM> centerA;

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

        const data_t dxA = x - (data_t)m_params.centerA[0];
        const data_t dyA = y - (data_t)m_params.centerA[1];
        const data_t dzA = z - (data_t)m_params.centerA[2];
        const data_t rA2 = dxA * dxA + dyA * dyA + dzA * dzA;

        data_t chi = 1.0;
        data_t h11 = 1.0, h12 = 0.0, h13 = 0.0, h22 = 1.0, h23 = 0.0, h33 = 1.0;

        const data_t r_sq_4 = 4.0 * rA2;
        const data_t den    = r_sq_4 + (data_t)b0_sq;
        const data_t frac   = r_sq_4 / den;
        chi = frac * frac;

        if (chi < (data_t)1.0e-10) chi = (data_t)1.0e-10;

        data_t phi = 0.0;
        data_t Pi = 0.0;

        const data_t eps2 = (data_t)1.0e-24;
        const data_t rA2_reg = simd_max(rA2, eps2);
        const data_t rA = sqrt(rA2_reg);
        const data_t argA = (rA - (data_t)b0_sq / (4.0 * rA)) / (data_t)b0;
        phi = (data_t)(1.0 / sqrt(4.0 * M_PI)) * atan(argA);

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
            lapse = 1.0 - (data_t)3.0 * log(chi);
        }
        else if (m_params.initial_lapse_type == 3)
        {
            lapse = chi;
        }
        else if (m_params.initial_lapse_type == 4)
        {
            const data_t core_radius = (data_t)(0.3 * b0);
            const data_t scaled_r = rA / core_radius;
            const data_t scaled_r2 = scaled_r * scaled_r;
            const data_t scaled_r4 = scaled_r2 * scaled_r2;
            const data_t scaled_r8 = scaled_r4 * scaled_r4;
            lapse = 1.0 - exp(-scaled_r8);
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