/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(BOOSTEDBH_HPP_)
#error "This file should only be included through BoostedBH.hpp"
#endif

#ifndef BOOSTEDBH_IMPL_HPP_
#define BOOSTEDBH_IMPL_HPP_

#include "BoostedBH.hpp"
#include "DimensionDefinitions.hpp"
#include <cmath>

inline BoostedBH::BoostedBH(params_t a_params) : m_params(a_params) {}

template <class data_t>
AMREX_GPU_DEVICE data_t
BoostedBH::psi_minus_one(Coordinates<data_t> coords) const
{
    const data_t r         = center_dist(coords);
    const data_t cos_theta = (coords.z - m_params.center[2]) / r;
    const data_t P_squared = std::pow(m_params.momentum[0], 2) +
                             std::pow(m_params.momentum[1], 2) +
                             std::pow(m_params.momentum[2], 2);
    return psi0(r) +
           P_squared * psi2(r, cos_theta) / (m_params.mass * m_params.mass);
}

template <class data_t>
AMREX_GPU_DEVICE Tensor<2, data_t>
BoostedBH::Aij(Coordinates<data_t> a_coords) const
{
    const data_t r = center_dist(a_coords);
    const Tensor<1, data_t> l{(a_coords.x - m_params.center[0]) / r,
                              (a_coords.y - m_params.center[1]) / r,
                              (a_coords.z - m_params.center[2]) / r};
    const data_t l_dot_p = l[0] * m_params.momentum[0] +
                           l[1] * m_params.momentum[1] +
                           l[2] * m_params.momentum[2];

    Tensor<2, data_t> out;

    FOR (i, j)
    {
        const double delta = (i == j) ? 1 : 0;
        out[i][j]          = 1.5 *
                    (m_params.momentum[i] * l[j] + m_params.momentum[j] * l[i] -
                     (delta - l[i] * l[j]) * l_dot_p) /
                    (r * r);
    }
    return out;
}

/* PRIVATE */

template <class data_t>
AMREX_GPU_DEVICE data_t
BoostedBH::center_dist(Coordinates<data_t> a_coords) const
{
    data_t r = std::sqrt(std::pow(a_coords.x - m_params.center[0], 2) +
                         std::pow(a_coords.y - m_params.center[1], 2) +
                         std::pow(a_coords.z - m_params.center[2], 2));

    double minimum_r = 1e-6;
    return simd_max(r, minimum_r);
}

template <class data_t>
AMREX_GPU_DEVICE data_t BoostedBH::psi0(data_t a_r) const
{
    return m_params.mass / (2 * a_r);
}

template <class data_t>
AMREX_GPU_DEVICE data_t BoostedBH::psi2(data_t a_r, data_t a_cos_theta) const
{
    return psi2_0(a_r) + psi2_2(a_r) * (1.5 * a_cos_theta * a_cos_theta - 0.5);
}

template <class data_t>
AMREX_GPU_DEVICE data_t BoostedBH::psi2_0(data_t a_r) const
{
    const data_t psi0_here    = psi0(a_r);
    const data_t psi0_sq_here = psi0_here * psi0_here;
    return std::pow(1 + psi0_here, -5) * (psi0_here / 8) *
           (psi0_sq_here * psi0_sq_here + 5 * psi0_here * psi0_sq_here +
            10 * psi0_sq_here + 10 * psi0_here + 5);
}

template <class data_t>
AMREX_GPU_DEVICE data_t BoostedBH::psi2_2(data_t a_r) const
{
    const data_t psi0_here    = psi0(a_r);
    const data_t psi0_sq_here = psi0_here * psi0_here;
    return 0.05 * std::pow(1 + psi0_here, -5) * psi0_sq_here *
               (84 * psi0_here * psi0_sq_here * psi0_sq_here +
                378 * psi0_sq_here * psi0_sq_here +
                658 * psi0_here * psi0_sq_here + 539 * psi0_sq_here +
                192 * psi0_here + 15) +
           4.2 * psi0_here * psi0_sq_here * log(psi0_here / (1 + psi0_here));
}

#endif /* BOOSTEDBH_IMPL_HPP_ */
