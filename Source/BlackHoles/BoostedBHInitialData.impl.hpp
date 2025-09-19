/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(BOOSTEDBHINITIALDATA_HPP_)
#error "This file should only be included through BoostedBHInitialData.hpp"
#endif

#ifndef BOOSTEDBHINITIALDATA_IMPL_HPP_
#define BOOSTEDBHINITIALDATA_IMPL_HPP_

#include "AlwaysInline.hpp"
#include "DimensionDefinitions.hpp"
#include <cmath>

AMREX_FORCE_INLINE BoostedBHInitialData::BoostedBHInitialData(params_t a_params)
    : m_params(a_params)
{
}

[[nodiscard]] AMREX_FORCE_INLINE AMREX_GPU_DEVICE amrex::Real
BoostedBHInitialData::psi_minus_one(Coordinates coords) const
{
    const amrex::Real r         = center_dist(coords);
    const amrex::Real cos_theta = (coords.z - m_params.center[2]) / r;
    const amrex::Real P_squared = std::pow(m_params.momentum[0], 2) +
                                  std::pow(m_params.momentum[1], 2) +
                                  std::pow(m_params.momentum[2], 2);
    return psi0(r) +
           P_squared * psi2(r, cos_theta) / (m_params.mass * m_params.mass);
}

[[nodiscard]] AMREX_FORCE_INLINE
    AMREX_GPU_DEVICE amrex::Array2D<amrex::Real, 0, 3, 0, 3>
    BoostedBHInitialData::Aij(Coordinates a_coords) const
{
    const amrex::Real r = center_dist(a_coords);
    const amrex::Array1D<amrex::Real, 0, 3> l{
        (a_coords.x - m_params.center[0]) / r,
        (a_coords.y - m_params.center[1]) / r,
        (a_coords.z - m_params.center[2]) / r};
    const amrex::Real l_dot_p = l(0) * m_params.momentum[0] +
                                l(1) * m_params.momentum[1] +
                                l(2) * m_params.momentum[2];

    amrex::Array2D<amrex::Real, 0, 3, 0, 3> out;

    FOR (i, j)
    {
        const double delta = (i == j) ? 1 : 0;
        out(i, j)          = 1.5 *
                    (m_params.momentum[i] * l(j) + m_params.momentum[j] * l(i) -
                     (delta - l(i) * l(j)) * l_dot_p) /
                    (r * r);
    }
    return out;
}

/* PRIVATE */

[[nodiscard]] AMREX_FORCE_INLINE AMREX_GPU_DEVICE amrex::Real
BoostedBHInitialData::center_dist(Coordinates a_coords) const
{
    amrex::Real r = std::sqrt(std::pow(a_coords.x - m_params.center[0], 2) +
                              std::pow(a_coords.y - m_params.center[1], 2) +
                              std::pow(a_coords.z - m_params.center[2], 2));

    return std::max(r, 1e-6);
}

[[nodiscard]] AMREX_FORCE_INLINE AMREX_GPU_DEVICE amrex::Real
BoostedBHInitialData::psi0(amrex::Real a_r) const
{
    return m_params.mass / (2 * a_r);
}

[[nodiscard]] AMREX_FORCE_INLINE AMREX_GPU_DEVICE amrex::Real
BoostedBHInitialData::psi2(amrex::Real a_r, amrex::Real a_cos_theta) const
{
    return psi2_0(a_r) + psi2_2(a_r) * (1.5 * a_cos_theta * a_cos_theta - 0.5);
}

[[nodiscard]] AMREX_FORCE_INLINE AMREX_GPU_DEVICE amrex::Real
BoostedBHInitialData::psi2_0(amrex::Real a_r) const
{
    const amrex::Real psi0_here    = psi0(a_r);
    const amrex::Real psi0_sq_here = psi0_here * psi0_here;
    return std::pow(1 + psi0_here, -5) * (psi0_here / 8) *
           (psi0_sq_here * psi0_sq_here + 5 * psi0_here * psi0_sq_here +
            10 * psi0_sq_here + 10 * psi0_here + 5);
}

[[nodiscard]] AMREX_FORCE_INLINE AMREX_GPU_DEVICE amrex::Real
BoostedBHInitialData::psi2_2(amrex::Real a_r) const
{
    const amrex::Real psi0_here    = psi0(a_r);
    const amrex::Real psi0_sq_here = psi0_here * psi0_here;
    return 0.05 * std::pow(1 + psi0_here, -5) * psi0_sq_here *
               (84 * psi0_here * psi0_sq_here * psi0_sq_here +
                378 * psi0_sq_here * psi0_sq_here +
                658 * psi0_here * psi0_sq_here + 539 * psi0_sq_here +
                192 * psi0_here + 15) +
           4.2 * psi0_here * psi0_sq_here * log(psi0_here / (1 + psi0_here));
}

#endif /* BOOSTEDBHINITIALDATA_IMPL_HPP_ */
