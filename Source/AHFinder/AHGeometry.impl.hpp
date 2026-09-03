/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(AHGEOMETRY_HPP_)
#error "This file should only be included through AHGeometry.hpp"
#endif

#ifndef AHGEOMETRY_IMPL_HPP_
#define AHGEOMETRY_IMPL_HPP_

#include "TensorAlgebra.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

// AHGeometry is not a template, so every out-of-class definition below must
// have inline linkage: this file is pulled in by AHGeometry.hpp, and without
// it a second translation unit including that header would give duplicate
// symbols at link time.

inline int AHGeometry::choose_n_rings(int num_particles)
{
    // Aim for ~twice as many longitudes as latitudes
    int n_rings = 1;
    int target  = std::max(1, (int)std::round(std::sqrt(num_particles / 2.0)));
    for (int n = 1; n <= num_particles; ++n)
    {
        if (num_particles % n == 0 &&
            std::abs(n - target) < std::abs(n_rings - target))
        {
            n_rings = n;
        }
    }

    return n_rings;
}

inline AHGeometry::AHGeometry(int num_particles,
                              const std::array<double, AMREX_SPACEDIM> &center,
                              double guess_radius)
    : m_num_particles(num_particles), m_n_rings(choose_n_rings(num_particles)),
      m_ring_size(num_particles / m_n_rings), m_center(center),
      m_guess_radius(guess_radius), m_dir_x(num_particles),
      m_dir_y(num_particles), m_dir_z(num_particles), m_dhdx(num_particles),
      m_dhdy(num_particles), m_dhdz(num_particles), m_d2h_xx(num_particles),
      m_d2h_yy(num_particles), m_d2h_zz(num_particles), m_d2h_xy(num_particles),
      m_d2h_xz(num_particles), m_d2h_yz(num_particles),
      m_dh_dtheta(num_particles), m_dh_dphi(num_particles)
{
    set_directions();
}

// Lay the grid points out in rings on a sphere, one unit direction each
inline void AHGeometry::set_directions()
{
    for (int i = 0; i < m_n_rings; ++i)
    {
        // theta() offsets by half a cell so we don't get points on the poles
        const double sin_theta = std::sin(theta(i));
        const double cos_theta = std::cos(theta(i));

        for (int j = 0; j < m_ring_size; ++j)
        {
            const double phi_ij = phi(j);
            const int idx       = i * m_ring_size + j;

            m_dir_x[idx] = std::cos(phi_ij) * sin_theta;
            m_dir_y[idx] = std::sin(phi_ij) * sin_theta;
            m_dir_z[idx] = cos_theta;
        }
    }
}

inline amrex::GpuArray<int, 4> AHGeometry::neighbours(int idx) const
{
    // Get the 4 neighbours of a particle (north, south, east and west).
    // If we are at the top/bottom ring, we can take the north/south
    // neighbour respectively as the opposite point on the same ring
    int ring_num = idx / m_ring_size;
    int ring_pos = idx % m_ring_size;

    int opposite_point =
        ring_num * m_ring_size + (ring_pos + m_ring_size / 2) % m_ring_size;

    int north, south, east, west;

    // North/south neighbours are next ring above/below.
    if (ring_num > 0)
    {
        north = (ring_num - 1) * m_ring_size + ring_pos;
    }
    else
    {
        north = opposite_point;
    }

    if (ring_num < m_n_rings - 1)
    {
        south = (ring_num + 1) * m_ring_size + ring_pos;
    }
    else
    {
        south = opposite_point;
    }

    // East/west neighbours are adjacent in the same ring
    // Modulo to wrap around
    east = ring_num * m_ring_size + (ring_pos + 1) % m_ring_size;
    west = ring_num * m_ring_size + (ring_pos - 1 + m_ring_size) % m_ring_size;

    return {north, south, east, west};
}

inline Tensor::Rank1 AHGeometry::direction(int idx) const
{
    return {m_dir_x[idx], m_dir_y[idx], m_dir_z[idx]};
}

AMREX_FORCE_INLINE amrex::GpuArray<amrex::Real, 3>
AHGeometry::ring_gradient(const double *field_ptr, int idx, double r) const
{
    // Perform finite differences with a particle and its
    // neighbours in spherical (theta/phi) coordinates and convert
    // them to cartesian (x/y/z)
    const amrex::GpuArray<int, 4> nb = neighbours(idx);

    amrex::Real d_dtheta =
        (field_ptr[nb[1]] - field_ptr[nb[0]]) / (2.0 * d_theta());
    amrex::Real d_dphi =
        (field_ptr[nb[2]] - field_ptr[nb[3]]) / (2.0 * d_phi());

    // Convert to Cartesian via the spherical Jacobian
    double sin_theta = std::sin(theta(idx / m_ring_size));
    double cos_theta = std::cos(theta(idx / m_ring_size));
    double sin_phi   = std::sin(phi(idx % m_ring_size));
    double cos_phi   = std::cos(phi(idx % m_ring_size));

    amrex::Real dx =
        d_dtheta * cos_theta * cos_phi / r - d_dphi * sin_phi / (r * sin_theta);
    amrex::Real dy =
        d_dtheta * cos_theta * sin_phi / r + d_dphi * cos_phi / (r * sin_theta);
    amrex::Real dz = -d_dtheta * sin_theta / r;

    return {dx, dy, dz};
}

inline void AHGeometry::set_h_derivatives(const std::vector<double> &h)
{
    h_derivs(h);
    h_hessian(h);
}

inline void AHGeometry::h_derivs(const std::vector<double> &h)
{
    // Calculate dh/dx, dh/dy, dh/dz for all particles.
    const double *h_ptr = h.data();
    double *dhdx_ptr    = m_dhdx.data();
    double *dhdy_ptr    = m_dhdy.data();
    double *dhdz_ptr    = m_dhdz.data();

    for (int id = 0; id < m_num_particles; ++id)
    {
        const amrex::GpuArray<amrex::Real, 3> grad_h =
            ring_gradient(h_ptr, id, h_ptr[id]);

        dhdx_ptr[id] = grad_h[0];
        dhdy_ptr[id] = grad_h[1];
        dhdz_ptr[id] = grad_h[2];
    }
}

inline void AHGeometry::h_hessian(const std::vector<double> &h)
{
    // Calculate the Hessian of h
    // Must be called after h_derivs() has populated m_dhdx/m_dhdy/m_dhdz.
    const double *h_ptr    = h.data();
    const double *dhdx_ptr = m_dhdx.data();
    const double *dhdy_ptr = m_dhdy.data();
    const double *dhdz_ptr = m_dhdz.data();

    double *d2h_xx_ptr = m_d2h_xx.data();
    double *d2h_yy_ptr = m_d2h_yy.data();
    double *d2h_zz_ptr = m_d2h_zz.data();
    double *d2h_xy_ptr = m_d2h_xy.data();
    double *d2h_xz_ptr = m_d2h_xz.data();
    double *d2h_yz_ptr = m_d2h_yz.data();

    for (int id = 0; id < m_num_particles; ++id)
    {
        const double r = h_ptr[id];

        // d/dx_j(dh/dx), d/dx_j(dh/dy), d/dx_j(dh/dz)
        const amrex::GpuArray<amrex::Real, 3> grad_x =
            ring_gradient(dhdx_ptr, id, r);
        const amrex::GpuArray<amrex::Real, 3> grad_y =
            ring_gradient(dhdy_ptr, id, r);
        const amrex::GpuArray<amrex::Real, 3> grad_z =
            ring_gradient(dhdz_ptr, id, r);

        d2h_xx_ptr[id] = grad_x[0];
        d2h_yy_ptr[id] = grad_y[1];
        d2h_zz_ptr[id] = grad_z[2];

        d2h_xy_ptr[id] = 0.5 * (grad_x[1] + grad_y[0]);
        d2h_xz_ptr[id] = 0.5 * (grad_x[2] + grad_z[0]);
        d2h_yz_ptr[id] = 0.5 * (grad_y[2] + grad_z[1]);
    }
}

inline Tensor::Rank1 AHGeometry::grad_h(int idx) const
{
    return {m_dhdx[idx], m_dhdy[idx], m_dhdz[idx]};
}

inline Tensor::Rank2 AHGeometry::hess_h(int idx) const
{
    Tensor::Rank2 hess_h_LL;
    hess_h_LL(0, 0) = m_d2h_xx[idx];
    hess_h_LL(1, 1) = m_d2h_yy[idx];
    hess_h_LL(2, 2) = m_d2h_zz[idx];
    hess_h_LL(0, 1) = hess_h_LL(1, 0) = m_d2h_xy[idx];
    hess_h_LL(0, 2) = hess_h_LL(2, 0) = m_d2h_xz[idx];
    hess_h_LL(1, 2) = hess_h_LL(2, 1) = m_d2h_yz[idx];

    return hess_h_LL;
}

// Closest pair of grid points anywhere on the surface: the rings pinch
// together towards the poles, so the phi spacing carries the sin(theta)
inline double AHGeometry::min_ring_spacing(const std::vector<double> &h) const
{
    double min_spacing = std::numeric_limits<double>::max();

    for (int i = 0; i < m_n_rings; ++i)
    {
        const double sin_theta = std::sin(theta(i));

        for (int j = 0; j < m_ring_size; ++j)
        {
            const double r = h[i * m_ring_size + j];

            min_spacing = std::min(min_spacing, r * d_theta());
            min_spacing = std::min(min_spacing, r * sin_theta * d_phi());
        }
    }

    return min_spacing;
}

inline void
AHGeometry::set_surface_data(const std::vector<double> *h,
                             const std::vector<Tensor::Rank2> *gamma_LL,
                             const std::vector<Tensor::Rank2> *K_LL,
                             const std::vector<double> *K)
{
    m_h        = h;
    m_gamma_LL = gamma_LL;
    m_K_LL     = K_LL;
    m_K        = K;
}

inline Tensor::Rank1 AHGeometry::s_L(int idx) const
{
    const Tensor::Rank1 n_L      = direction(idx);
    const Tensor::Rank1 grad_h_L = grad_h(idx);

    Tensor::Rank1 F_L;
    FOR (i)
    {
        F_L(i) = n_L(i) - grad_h_L(i);
    }

    const Tensor::Rank2 gamma_UU =
        TensorAlgebra::compute_inverse((*m_gamma_LL)[idx]);
    const amrex::Real lambda =
        std::sqrt(TensorAlgebra::compute_dot_product(F_L, F_L, gamma_UU));

    Tensor::Rank1 sL;
    FOR (i)
    {
        sL(i) = F_L(i) / lambda;
    }

    return sL;
}

inline Tensor::Rank1 AHGeometry::s_U(int idx) const
{
    const Tensor::Rank2 gamma_UU =
        TensorAlgebra::compute_inverse((*m_gamma_LL)[idx]);

    return TensorAlgebra::raise_all(s_L(idx), gamma_UU);
}

inline Tensor::Rank2 AHGeometry::pi_LL(int idx) const
{
    const Tensor::Rank2 &K_LL_idx     = (*m_K_LL)[idx];
    const Tensor::Rank2 &gamma_LL_idx = (*m_gamma_LL)[idx];
    const amrex::Real K_idx           = (*m_K)[idx];

    Tensor::Rank2 pi;
    FOR (i, j)
    {
        pi(i, j) = K_LL_idx(i, j) - K_idx * gamma_LL_idx(i, j);
    }

    return pi;
}

// h(theta, phi) is already a scalar field native to the ring grid
// dh/dtheta, dh/dphi directly calculated
inline void AHGeometry::set_h_gradients()
{
    for (int i = 0; i < m_n_rings; ++i)
    {
        for (int j = 0; j < m_ring_size; ++j)
        {
            const int ij = i * m_ring_size + j;

            const amrex::GpuArray<int, 4> nb = neighbours(ij);
            const int north                  = nb[0];
            const int south                  = nb[1];
            const int east                   = nb[2];
            const int west                   = nb[3];

            m_dh_dtheta[ij] =
                ((*m_h)[south] - (*m_h)[north]) / (2.0 * d_theta());
            m_dh_dphi[ij] = ((*m_h)[east] - (*m_h)[west]) / (2.0 * d_phi());
        }
    }
}

// d x^i/d theta = e_theta^i
// where x^i are cartesian coordiantes of the surface
inline Tensor::Rank1 AHGeometry::e_theta(int i, int j) const
{
    const int ij = i * m_ring_size + j;

    const amrex::Real cos_phi   = std::cos(phi(j));
    const amrex::Real sin_phi   = std::sin(phi(j));
    const amrex::Real cos_theta = std::cos(theta(i));
    const amrex::Real sin_theta = std::sin(theta(i));

    // e_theta^i = d(x^i)/dtheta
    Tensor::Rank1 dx_dth;
    dx_dth(0) = m_dh_dtheta[ij] * cos_phi * sin_theta +
                (*m_h)[ij] * cos_phi * cos_theta;
    dx_dth(1) = m_dh_dtheta[ij] * sin_phi * sin_theta +
                (*m_h)[ij] * sin_phi * cos_theta;
    dx_dth(2) = m_dh_dtheta[ij] * cos_theta - (*m_h)[ij] * sin_theta;

    return dx_dth;
}

// d x^i/d phi = e_phi^i
// where x^i are cartesian coordiantes of the surface
inline Tensor::Rank1 AHGeometry::e_phi(int i, int j) const
{
    const int ij = i * m_ring_size + j;

    const amrex::Real cos_phi   = std::cos(phi(j));
    const amrex::Real sin_phi   = std::sin(phi(j));
    const amrex::Real cos_theta = std::cos(theta(i));
    const amrex::Real sin_theta = std::sin(theta(i));

    // e_phi^i = d(x^i)/dphi
    Tensor::Rank1 dx_dph;
    dx_dph(0) =
        m_dh_dphi[ij] * cos_phi * sin_theta - (*m_h)[ij] * sin_phi * sin_theta;
    dx_dph(1) =
        m_dh_dphi[ij] * sin_phi * sin_theta + (*m_h)[ij] * cos_phi * sin_theta;
    dx_dph(2) = m_dh_dphi[ij] * cos_theta;

    return dx_dph;
}

// q_ab is the 2x2 metric adapted to the surface
// calculated from the pullback onto the surface of the 3 metric
// the indices ab are [theta,phi]
inline QAB AHGeometry::q_ab(int i, int j) const
{
    const int ij = i * m_ring_size + j;

    const Tensor::Rank2 &gamma = (*m_gamma_LL)[ij];

    const Tensor::Rank1 e_th = e_theta(i, j);
    const Tensor::Rank1 e_ph = e_phi(i, j);

    // q_ab = gamma_ij e_a^i e_b^j
    QAB qLL;
    qLL(0, 0) = TensorAlgebra::compute_dot_product(e_th, e_th, gamma);
    qLL(0, 1) = TensorAlgebra::compute_dot_product(e_th, e_ph, gamma);
    qLL(1, 0) = qLL(0, 1);
    qLL(1, 1) = TensorAlgebra::compute_dot_product(e_ph, e_ph, gamma);

    return qLL;
}

// area scale factor of surface
// root_q is the square root of det q_ab
inline amrex::Real AHGeometry::root_det_q(int i, int j) const
{
    const QAB qLL = q_ab(i, j);

    const amrex::Real det_q = qLL(0, 0) * qLL(1, 1) - qLL(0, 1) * qLL(1, 0);

    return std::sqrt(det_q);
}

// integrtes : Area = integral dA
// dA = root_q * dtheta * dphi
inline amrex::Real AHGeometry::area()
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_h != nullptr && m_gamma_LL != nullptr,
        "AHGeometry::area() called before set_surface_data()");

    // Refresh dh/dtheta, dh/dphi from the currently bound h so the tangent
    // vectors below can't be built from a stale surface
    set_h_gradients();

    amrex::Real total = 0.0;
    for (int i = 0; i < m_n_rings; ++i)
    {
        for (int j = 0; j < m_ring_size; ++j)
        {
            total += root_det_q(i, j);
        }
    }

    return total * d_theta() * d_phi();
}

// P_k = (1/8pi) * integral pi_kj s^j dA, for each Cartesian direction k
inline Tensor::Rank1 AHGeometry::momentum()
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_h != nullptr && m_gamma_LL != nullptr && m_K_LL != nullptr &&
            m_K != nullptr,
        "AHGeometry::momentum() called before set_surface_data()");

    // Refresh dh/dtheta, dh/dphi from the currently bound h so the tangent
    // vectors used by q_ab() below can't be built from a stale surface
    set_h_gradients();

    Tensor::Rank1 total;
    FOR (k)
    {
        total(k) = 0.0;
    }

    for (int i = 0; i < m_n_rings; ++i)
    {
        for (int j = 0; j < m_ring_size; ++j)
        {
            const int idx = i * m_ring_size + j;

            const Tensor::Rank2 pi   = pi_LL(idx);
            const Tensor::Rank1 sU   = s_U(idx);
            const amrex::Real root_q = root_det_q(i, j);

            FOR (k)
            {
                amrex::Real pi_k_dot_s = 0.0;
                FOR (l)
                {
                    pi_k_dot_s += pi(k, l) * sU(l);
                }
                total(k) += pi_k_dot_s * root_q;
            }
        }
    }

    // the 1/8 pi here is from the komar integral for momentum
    FOR (k)
    {
        total(k) *= d_theta() * d_phi() / (8.0 * M_PI);
    }

    return total;
}

#endif /* AHGEOMETRY_IMPL_HPP_ */
