/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#include "AHGeometry.hpp"

#include "TensorAlgebra.hpp"

#include <cmath>

void AHGeometry::refresh(int num_theta, int num_phi, const std::vector<double> *f,
                         const std::vector<Tensor::Rank2> *gamma_LL)
{
    m_num_theta = num_theta;
    m_num_phi   = num_phi;
    m_f         = f;
    m_gamma_LL  = gamma_LL;

    const amrex::Real d_theta = M_PI / m_num_theta;
    const amrex::Real d_phi   = 2.0 * M_PI / m_num_phi;
    m_d_theta_d_phi           = d_theta * d_phi;

    set_f_gradients();
}


// used for finite difference derivatives on the surface
// assigns 'off the grid' points for derivatives at the edge
std::array<int, 4> AHGeometry::neighbours(int i, int j) const
{
    const int opposite_point = i * m_num_phi + (j + m_num_phi / 2) % m_num_phi;

    int north=0; 
    int south=0; 
    int east=0; 
    int west=0;

    if (i > 0)
        north = (i - 1) * m_num_phi + j;
    else
        north = opposite_point;

    if (i < m_num_theta - 1)
        south = (i + 1) * m_num_phi + j;
    else
        south = opposite_point;

    east = i * m_num_phi + (j + 1) % m_num_phi;
    west = i * m_num_phi + (j - 1 + m_num_phi) % m_num_phi;

    return {north, south, east, west};
}

// f(theta, phi) is already a scalar field native to the ring grid
// df/dtheta, df/dphi directly calcaulted
void AHGeometry::set_f_gradients()
{
    const int num_particles = m_num_theta * m_num_phi;
    m_df_dtheta.resize(num_particles);
    m_df_dphi.resize(num_particles);

    const amrex::Real d_theta = M_PI / m_num_theta;
    const amrex::Real d_phi   = 2.0 * M_PI / m_num_phi;

    for (int i = 0; i < m_num_theta; ++i)
    {
        for (int j = 0; j < m_num_phi; ++j)
        {
            const int ij = i * m_num_phi + j;

            const std::array<int, 4> nb = neighbours(i, j);
            const int north = nb[0];
            const int south = nb[1];
            const int east  = nb[2];
            const int west  = nb[3];

            m_df_dtheta[ij] = ((*m_f)[south] - (*m_f)[north]) / (2.0 * d_theta);
            m_df_dphi[ij]   = ((*m_f)[east] - (*m_f)[west]) / (2.0 * d_phi);
        }
    }
}


// d x^i/d theta = e_theta^i
// where x^i are cartesian coordiantes of the surface
Tensor::Rank1 AHGeometry::e_theta(int i, int j) const
{
    const int ij = i * m_num_phi + j;

    const amrex::Real theta = (i + 0.5) * M_PI / m_num_theta;
    const amrex::Real phi   = j * 2.0 * M_PI / m_num_phi;

    const amrex::Real cos_phi   = std::cos(phi);
    const amrex::Real sin_phi   = std::sin(phi);
    const amrex::Real cos_theta = std::cos(theta);
    const amrex::Real sin_theta = std::sin(theta);

    // e_theta^i = d(x^i)/dtheta
    Tensor::Rank1 dx_dth;
    dx_dth(0) = m_df_dtheta[ij] * cos_phi * sin_theta
         + (*m_f)[ij] * cos_phi * cos_theta;
    dx_dth(1) = m_df_dtheta[ij] * sin_phi * sin_theta
         + (*m_f)[ij] * sin_phi * cos_theta;
    dx_dth(2) = m_df_dtheta[ij] * cos_theta
         - (*m_f)[ij] * sin_theta;

    return dx_dth;
}

// d x^i/d phi = e_phi^i
// where x^i are cartesian coordiantes of the surface
Tensor::Rank1 AHGeometry::e_phi(int i, int j) const
{
    const int ij = i * m_num_phi + j;

    const amrex::Real theta = (i + 0.5) * M_PI / m_num_theta;
    const amrex::Real phi   = j * 2.0 * M_PI / m_num_phi;

    const amrex::Real cos_phi   = std::cos(phi);
    const amrex::Real sin_phi   = std::sin(phi);
    const amrex::Real cos_theta = std::cos(theta);
    const amrex::Real sin_theta = std::sin(theta);

    // e_phi^i = d(x^i)/dphi
    Tensor::Rank1 dx_dph;
    dx_dph(0) = m_df_dphi[ij] * cos_phi * sin_theta
         - (*m_f)[ij] * sin_phi * sin_theta;
    dx_dph(1) = m_df_dphi[ij] * sin_phi * sin_theta
         + (*m_f)[ij] * cos_phi * sin_theta;
    dx_dph(2) = m_df_dphi[ij] * cos_theta;

    return dx_dph;
}

// q_ab is the 2x2 metric adapted to the surface
// calculated from the pullback onto the surface of the 3 metric
// the indices ab are [theta,phi]
QAB AHGeometry::q_ab(int i, int j) const
{
    const int ij = i * m_num_phi + j;

    const Tensor::Rank2 &gamma = (*m_gamma_LL)[ij];

    const Tensor::Rank1 e_th = e_theta(i, j);
    const Tensor::Rank1 e_ph = e_phi(i, j);

    // q_ab = gamma_ij e_a^i e_b^j
    QAB qLL;
    qLL(0, 0) = TensorAlgebra::compute_dot_product(e_th, e_th, gamma);
    qLL(0, 1) = TensorAlgebra::compute_dot_product(e_th, e_ph, gamma);
    qLL(1, 0) = q(0, 1);
    qLL(1, 1) = TensorAlgebra::compute_dot_product(e_ph, e_ph, gamma);

    return qLL;
}

// area scale factor of surface
// root_q is the square root of det q_ab
amrex::Real AHGeometry::root_det_q(int i, int j) const
{
    const QAB q = q_ab(i, j);

    const amrex::Real det_q = q(0, 0) * q(1, 1) - q(0, 1) * q(1, 0);

    return std::sqrt(det_q);
}


// integrtes : Area = integral dA
// dA = root_q * dtheta * dphi
amrex::Real AHGeometry::area() const
{
    amrex::Real total = 0.0;
    for (int i = 0; i < m_num_theta; ++i)
    {
        for (int j = 0; j < m_num_phi; ++j)
        {
            total += root_det_q(i, j);
        }
    }

    return total * m_d_theta_d_phi;
}
