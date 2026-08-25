/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(OSCILLATONINITIALDATA_HPP_)
#error "This file should only be included through OscillatonInitialData.hpp"
#endif

#ifndef OSCILLATONINITIALDATA_IMPL_HPP_
#define OSCILLATONINITIALDATA_IMPL_HPP_

#include "OscillatonInitialData.hpp"
#include "TensorAlgebra.hpp"

#include <cmath>

AMREX_FORCE_INLINE
OscillatonInitialData::OscillatonInitialData(params_t a_params,
                                             amrex::Real a_dx)
    : m_dx(a_dx), m_params(a_params)
{
}

template <std::size_t N>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
OscillatonInitialData::evaluate_chebyshev(
    amrex::Real x, const std::array<amrex::Real, N> &coefficients) const
{
    // Clenshaw recurrence for sum_k coefficients[k] T_k(x).
    amrex::Real b_k_plus_1 = 0.0;
    amrex::Real b_k_plus_2 = 0.0;
    for (int k = static_cast<int>(N) - 1; k >= 1; --k)
    {
        const amrex::Real b_k = 2.0 * x * b_k_plus_1 - b_k_plus_2 +
                                coefficients[static_cast<std::size_t>(k)];
        b_k_plus_2 = b_k_plus_1;
        b_k_plus_1 = b_k;
    }
    return coefficients[0] + x * b_k_plus_1 - b_k_plus_2;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
OscillatonInitialData::operator()(int ix, int iy, int iz,
                                  const amrex::Array4<amrex::Real> &state) const
{
    const amrex::CellData<amrex::Real> &state_cell_data =
        state.cellData(ix, iy, iz);
    const Coordinates coords(amrex::IntVect(ix, iy, iz), m_dx, m_params.center);

    const std::array<amrex::Real, GR_SPACEDIM> position{coords.x, coords.y,
                                                        coords.z};
    const amrex::Real r2 =
        coords.x * coords.x + coords.y * coords.y + coords.z * coords.z;

    // Compactified coordinate used by the metric and lapse fits.
    const amrex::Real geometry_scale2      = geometry_scale * geometry_scale;
    const amrex::Real geometry_denominator = r2 + geometry_scale2;
    const amrex::Real geometry_s           = r2 / geometry_denominator;
    const amrex::Real geometry_root =
        geometry_scale / std::sqrt(geometry_denominator);
    const amrex::Real geometry_x = 2.0 * geometry_s - 1.0;

    const amrex::Real compactness_polynomial =
        evaluate_chebyshev(geometry_x, m_compactness_coefficients);
    const amrex::Real compactness =
        geometry_s * geometry_root * compactness_polynomial;
    const amrex::Real grr = 1.0 / (1.0 - compactness);
    const amrex::Real chi = std::pow(grr, -1.0 / 3.0);

    const amrex::Real lapse_correction =
        geometry_root * geometry_root * geometry_root *
        evaluate_chebyshev(geometry_x, m_lapse_correction_coefficients);
    const amrex::Real lapse =
        std::sqrt(1.0 - compactness) * std::exp(lapse_correction);

    // The scalar profile uses a different compactification scale.
    const amrex::Real scalar_scale2      = scalar_scale * scalar_scale;
    const amrex::Real scalar_denominator = r2 + scalar_scale2;
    const amrex::Real scalar_s           = r2 / scalar_denominator;
    const amrex::Real scalar_root =
        scalar_scale / std::sqrt(scalar_denominator);
    const amrex::Real scalar_x = 2.0 * scalar_s - 1.0;
    const amrex::Real scalar_polynomial =
        evaluate_chebyshev(scalar_x, m_scalar_exponent_coefficients);
    const amrex::Real scalar_momentum =
        scalar_central_value *
        std::exp(-(scalar_s / scalar_root) * scalar_polynomial);

    state_cell_data[c_chi]   = chi;
    state_cell_data[c_lapse] = lapse;
    state_cell_data[c_Pi]    = scalar_momentum;

    // Evaluate (g_rr-1)/r^2 without a cancellation or a special case at r=0.
    const amrex::Real radial_factor =
        geometry_root * compactness_polynomial /
        (geometry_denominator * (1.0 - compactness));
    FOR2_SYM(i, j)
    {
        const amrex::Real gamma_ij = TensorAlgebra::delta(i, j) +
                                     radial_factor * position[i] * position[j];
        state_cell_data[sym_var_idx(c_h11, i, j)] = chi * gamma_ij;
    }

    // This will need to be set by the GammaCalculator, since the metric
    // is not conformally flat, but we will zero it here for now.
    FOR (i)
    {
        state_cell_data[c_Gamma1 + i] = 0.0;
    }

    // phi, K, A_ij, Theta, shift and B are all zero at this phase.
}

#endif /* OSCILLATONINITIALDATA_IMPL_HPP_ */
