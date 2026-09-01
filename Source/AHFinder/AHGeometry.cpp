/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#include "AHGeometry.hpp"

#include "TensorAlgebra.hpp"

#include <cmath>

void AHGeometry::refresh(int num_theta, int num_phi, const std::vector<double> *f,
                         const std::vector<Tensor::Rank2> *gamma_LL,
                         const std::vector<double> *dhdx,
                         const std::vector<double> *dhdy,
                         const std::vector<double> *dhdz)
{
    m_num_theta = num_theta;
    m_num_phi   = num_phi;
    m_f         = f;
    m_gamma_LL  = gamma_LL;
    m_dhdx      = dhdx;
    m_dhdy      = dhdy;
    m_dhdz      = dhdz;

    const amrex::Real d_theta = M_PI / m_num_theta;
    const amrex::Real d_phi   = 2.0 * M_PI / m_num_phi;
    m_d_theta_d_phi           = d_theta * d_phi;

    set_f_gradients();
}

void AHGeometry::set_f_gradients()
{
    const int num_particles = m_num_theta * m_num_phi;
    m_df_dtheta.resize(num_particles);
    m_df_dphi.resize(num_particles);

    for (int i = 0; i < m_num_theta; ++i)
    {
        for (int j = 0; j < m_num_phi; ++j)
        {
            const int ij = i * m_num_phi + j;

            const amrex::Real theta = (i + 0.5) * M_PI / m_num_theta;
            const amrex::Real phi   = j * 2.0 * M_PI / m_num_phi;

            const amrex::Real cos_theta = std::cos(theta);
            const amrex::Real sin_theta = std::sin(theta);
            const amrex::Real cos_phi   = std::cos(phi);
            const amrex::Real sin_phi   = std::sin(phi);

            const amrex::Real dhdx = (*m_dhdx)[ij];
            const amrex::Real dhdy = (*m_dhdy)[ij];
            const amrex::Real dhdz = (*m_dhdz)[ij];

            const amrex::Real f = (*m_f)[ij];

            // Project the Cartesian gradient of h onto the theta/phi
            // unit directions (inverse of ring_gradient's Jacobian).
            m_df_dtheta[ij] = f * (cos_theta * cos_phi * dhdx +
                                   cos_theta * sin_phi * dhdy -
                                   sin_theta * dhdz);
            m_df_dphi[ij]   = f * sin_theta * (-sin_phi * dhdx + cos_phi * dhdy);
        }
    }
}

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
    Tensor::Rank1 e;
    e(0) = m_df_dtheta[ij] * cos_phi * sin_theta
         + (*m_f)[ij] * cos_phi * cos_theta;
    e(1) = m_df_dtheta[ij] * sin_phi * sin_theta
         + (*m_f)[ij] * sin_phi * cos_theta;
    e(2) = m_df_dtheta[ij] * cos_theta
         - (*m_f)[ij] * sin_theta;

    return e;
}

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
    Tensor::Rank1 e;
    e(0) = m_df_dphi[ij] * cos_phi * sin_theta
         - (*m_f)[ij] * sin_phi * sin_theta;
    e(1) = m_df_dphi[ij] * sin_phi * sin_theta
         + (*m_f)[ij] * cos_phi * sin_theta;
    e(2) = m_df_dphi[ij] * cos_theta;

    return e;
}

QAB AHGeometry::q_ab(int i, int j) const
{
    const int ij = i * m_num_phi + j;

    const Tensor::Rank2 &gamma = (*m_gamma_LL)[ij];

    const Tensor::Rank1 e_th = e_theta(i, j);
    const Tensor::Rank1 e_ph = e_phi(i, j);

    // q_ab = gamma_ij e_a^i e_b^j
    QAB q;
    q(0, 0) = TensorAlgebra::compute_dot_product(e_th, e_th, gamma);
    q(0, 1) = TensorAlgebra::compute_dot_product(e_th, e_ph, gamma);
    q(1, 0) = q(0, 1);
    q(1, 1) = TensorAlgebra::compute_dot_product(e_ph, e_ph, gamma);

    return q;
}

amrex::Real AHGeometry::root_det_q(int i, int j) const
{
    const QAB q = q_ab(i, j);

    const amrex::Real det_q = q(0, 0) * q(1, 1) - q(0, 1) * q(1, 0);

    return std::sqrt(det_q);
}

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
