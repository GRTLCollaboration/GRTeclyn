/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef AHGEOMETRY_HPP_
#define AHGEOMETRY_HPP_

#include "Tensor.hpp"
#include <AMReX_REAL.H>
#include <array>
#include <vector>

// 2x2 tensor type for the surface's induced 2-metric q_ab (a, b index
// {theta, phi}); Tensor.hpp has no premade alias for a non-3-dimensional
// rank-2 tensor, so this is defined locally.
using QAB = Tensor::GeneralRank<2, 2, 2>;

// Computes geometric diagnostics of the AH surface (e.g. proper area).
// AHFinder points AHGeometry at its persistent per-particle arrays via
// refresh(), which should be called once AHFinder's surface data (f,
// gamma_ij) holds meaningful values -- e.g. after the pseudo-time
// evolution in AHFinder::find() has converged. AHGeometry derives
// df_dtheta/df_dphi from f itself; everything else is read live
// through the stored pointers, so no repeated refresh() calls are
// needed just to keep pace with AHFinder's per-step state.
class AHGeometry
{
  private:
    // Ring (latitude x longitude) grid dimensions, copied from AHFinder.
    int m_num_theta = 0;
    int m_num_phi   = 0;

    // Product of the (constant) ring-grid cell widths dtheta * dphi,
    // recomputed whenever the grid dimensions are refreshed.
    amrex::Real m_d_theta_d_phi = 0.0;

    // Derived surface data at each ring-grid point, flat-indexed as
    // i * m_num_phi + j to match AHFinder's ring-grid particle order.
    // Owned by AHGeometry: recomputed by set_f_gradients() each refresh.
    std::vector<amrex::Real> m_df_dtheta{};
    std::vector<amrex::Real> m_df_dphi{};

    // Pointers to AHFinder's persistent surface radius f (= h) and
    // physical 3-metric gamma_ij (same flat indexing as above). Not
    // owned by AHGeometry: set once by refresh() and read fresh each
    // time they're needed, so they always reflect AHFinder's latest
    // computed values.
    const std::vector<double> *m_f = nullptr;
    const std::vector<Tensor::Rank2> *m_gamma_LL = nullptr;

    // Flat indices of the 4 ring-grid neighbours of (i, j), ordered
    // {north, south, east, west}. Mirrors AHFinder::neighbours(): phi
    // wraps periodically, and at the pole-adjacent rings (where no
    // north/south ring exists) the missing neighbour is replaced by
    // the point at the same theta but phi + pi -- the exact reflection
    // across the pole, since a smooth field's (theta, phi) values
    // satisfy f(-theta, phi) = f(theta, phi + pi) there.
    std::array<int, 4> neighbours(int i, int j) const;

    // Fills m_df_dtheta, m_df_dphi at every ring-grid point via a
    // central difference of m_f over its ring-grid neighbours. f is
    // already a scalar field native to the (theta, phi) grid -- no
    // Cartesian round-trip needed to get its theta/phi derivatives.
    // Called from refresh().
    void set_f_gradients();

    // Tangent vector e_theta^i = d(x^i)/dtheta (Cartesian components)
    // at ring-grid point (i, j), from m_f and m_df_dtheta there. Used
    // only by q_ab().
    Tensor::Rank1 e_theta(int i, int j) const;

    // Tangent vector e_phi^i = d(x^i)/dphi (Cartesian components) at
    // ring-grid point (i, j), from m_f and m_df_dphi there. Used only
    // by q_ab().
    Tensor::Rank1 e_phi(int i, int j) const;

  public:
    AHGeometry() = default;

    // Copy over the current ring-grid dimensions, point at AHFinder's
    // persistent surface radius f and physical 3-metric, and derive
    // df_dtheta/df_dphi from f via set_f_gradients().
    void refresh(int num_theta, int num_phi, const std::vector<double> *f,
                 const std::vector<Tensor::Rank2> *gamma_LL);

    // Induced 2-metric q_ab of the surface at ring-grid point (i, j),
    // with a, b indexing {theta, phi}: q_ab = gamma_ij e_a^i e_b^j.
    QAB q_ab(int i, int j) const;

    // sqrt(det q_ab) at ring-grid point (i, j).
    amrex::Real root_det_q(int i, int j) const;

    // Total proper area of the surface: the sum of
    // root_det_q(i, j) * dtheta * dphi over the whole ring grid.
    amrex::Real area() const;
};

#endif /* AHGEOMETRY_HPP_ */
