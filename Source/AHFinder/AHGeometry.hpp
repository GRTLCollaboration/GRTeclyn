/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef AHGEOMETRY_HPP_
#define AHGEOMETRY_HPP_

#include "Tensor.hpp"
#include <AMReX_Array.H>
#include <AMReX_BLassert.H>
#include <AMReX_Extension.H>
#include <AMReX_REAL.H>
#include <array>
#include <cmath>
#include <vector>

// 2x2 tensor type for the surface's induced 2-metric q_ab (a, b index
// {theta, phi}); Tensor.hpp has no premade alias for a non-3-dimensional
// rank-2 tensor, so this is defined locally.
using QAB = Tensor::GeneralRank<2, 2, 2>;

// Owns the ring (latitude x longitude) grid the AH surface is discretised
// on, and everything derived from it: the fixed unit direction of each
// grid point, the finite-difference stencil (neighbours()), the Cartesian
// gradient and Hessian of the surface radius h, and the geometric
// diagnostics of the surface (induced 2-metric, proper area).
//
// The grid is built once, in the constructor, from the requested number of
// particles, so the directions and the stencil are usable immediately.
// AHFinder additionally binds AHGeometry to its persistent h and gamma_ij
// arrays with set_surface_data(); those are read live through the stored
// pointers by the area diagnostics, so they always see AHFinder's latest
// values without a copy. The evolution-facing entry points
// (set_h_derivatives(), min_ring_spacing()) instead take h by argument,
// because the AMReX time integrator evaluates them on temporary states
// that are not AHFinder's m_state.h.
class AHGeometry
{
  private:
    // Ring (latitude x longitude) decomposition of the particles, computed
    // from num_particles in the constructor.
    // m_n_rings * m_ring_size == m_num_particles.
    int m_num_particles;
    int m_n_rings;
    int m_ring_size;

    // Centre the surface is parameterised about, and the radius of the
    // initial (spherical) guess.
    std::array<double, AMREX_SPACEDIM> m_center;
    double m_guess_radius;

    // Fixed unit direction from the centre to each ring-grid point, flat
    // indexed as i * m_ring_size + j. Set once, in the constructor: the
    // surface point at idx is center() + h[idx] * direction(idx).
    std::vector<double> m_dir_x{};
    std::vector<double> m_dir_y{};
    std::vector<double> m_dir_z{};

    // Cartesian gradient and Hessian of the surface radius h, filled by
    // set_h_derivatives() and read back through grad_h() / hess_h().
    std::vector<double> m_dhdx{};
    std::vector<double> m_dhdy{};
    std::vector<double> m_dhdz{};
    std::vector<double> m_d2h_xx{};
    std::vector<double> m_d2h_yy{};
    std::vector<double> m_d2h_zz{};
    std::vector<double> m_d2h_xy{};
    std::vector<double> m_d2h_xz{};
    std::vector<double> m_d2h_yz{};

    // dh/dtheta and dh/dphi on the ring grid, refreshed by
    // set_h_gradients() at the start of area(). Distinct from m_dhdx etc.
    // above: these are the derivatives in the surface's own (theta, phi)
    // coordinates, and are used only by e_theta() / e_phi().
    std::vector<amrex::Real> m_dh_dtheta{};
    std::vector<amrex::Real> m_dh_dphi{};

    // Pointers to AHFinder's persistent surface radius h, physical
    // 3-metric gamma_ij, extrinsic curvature K_ij, and its trace K (same
    // flat indexing as above). Not owned by AHGeometry: set once by
    // set_surface_data() and read fresh each time they're needed, so they
    // always reflect AHFinder's latest values.
    const std::vector<double> *m_h               = nullptr;
    const std::vector<Tensor::Rank2> *m_gamma_LL = nullptr;
    const std::vector<Tensor::Rank2> *m_K_LL     = nullptr;
    const std::vector<double> *m_K               = nullptr;

    // Number of latitude rings: the divisor of num_particles closest to
    // sqrt(num_particles / 2), i.e. aiming for ~twice as many longitudes
    // as latitudes.
    static int choose_n_rings(int num_particles);

    // Fill m_dir_x/y/z with the unit vector at each ring-grid point.
    // Called from the constructor.
    void set_directions();

    // Differentiates a scalar field defined on the ring grid at flat index
    // idx, where the surface radius is r, and returns the result in
    // Cartesian components: central differences in (theta, phi) over
    // neighbours(idx), converted with the spherical Jacobian.
    [[nodiscard]] AMREX_FORCE_INLINE amrex::GpuArray<amrex::Real, 3>
    ring_gradient(const double *field_ptr, int idx, double r) const;

    // Fills m_dhdx/m_dhdy/m_dhdz from h. Called by set_h_derivatives().
    void h_derivs(const std::vector<double> &h);

    // Fills m_d2h_* by differentiating m_dhdx/m_dhdy/m_dhdz, so must run
    // after h_derivs(). Called by set_h_derivatives().
    void h_hessian(const std::vector<double> &h);

    // Fills m_dh_dtheta, m_dh_dphi at every ring-grid point via a central
    // difference of *m_h over its ring-grid neighbours. h is already a
    // scalar field native to the (theta, phi) grid -- no Cartesian
    // round-trip needed to get its theta/phi derivatives. Called from
    // area().
    void set_h_gradients();

    // Tangent vector e_theta^i = d(x^i)/dtheta (Cartesian components) at
    // ring-grid point (i, j), from *m_h and m_dh_dtheta there. Used only
    // by q_ab().
    [[nodiscard]] Tensor::Rank1 e_theta(int i, int j) const;

    // Tangent vector e_phi^i = d(x^i)/dphi (Cartesian components) at
    // ring-grid point (i, j), from *m_h and m_dh_dphi there. Used only by
    // q_ab().
    [[nodiscard]] Tensor::Rank1 e_phi(int i, int j) const;

  public:
    AHGeometry(int num_particles,
               const std::array<double, AMREX_SPACEDIM> &center,
               double guess_radius);

    // Ring grid

    [[nodiscard]] int num_particles() const { return m_num_particles; }
    [[nodiscard]] int n_rings() const { return m_n_rings; }
    [[nodiscard]] int ring_size() const { return m_ring_size; }

    // The (constant) ring-grid cell widths.
    [[nodiscard]] amrex::Real d_theta() const { return M_PI / m_n_rings; }
    [[nodiscard]] amrex::Real d_phi() const { return 2.0 * M_PI / m_ring_size; }

    // Polar angle of ring i, offset by half a cell so no grid point lands
    // on a pole, and azimuthal angle of position j within a ring.
    [[nodiscard]] amrex::Real theta(int i) const
    {
        return (i + 0.5) * d_theta();
    }
    [[nodiscard]] amrex::Real phi(int j) const { return j * d_phi(); }

    // Flat indices of the 4 ring-grid neighbours of the point at flat
    // index idx, ordered {north, south, east, west}. phi wraps
    // periodically, and at the pole-adjacent rings (where no north/south
    // ring exists) the missing neighbour is replaced by the point opposite
    // on the same ring
    [[nodiscard]] amrex::GpuArray<int, 4> neighbours(int idx) const;

    // Surface parameterisation

    [[nodiscard]] const std::array<double, AMREX_SPACEDIM> &center() const
    {
        return m_center;
    }

    [[nodiscard]] double guess_radius() const { return m_guess_radius; }

    // Unit direction from the centre to the ring-grid point at flat index
    // idx; the surface point there is center() + h[idx] * direction(idx).
    [[nodiscard]] Tensor::Rank1 direction(int idx) const;

    // Derivatives of the surface radius

    // Recomputes the Cartesian gradient and Hessian of h over the whole
    // ring grid. Must be called before grad_h() / hess_h() are read for a
    // given h; AHFinder::compute_theta() does so on every evaluation.
    void set_h_derivatives(const std::vector<double> &h);

    // (dh/dx, dh/dy, dh/dz) at flat index idx, from the last
    // set_h_derivatives().
    [[nodiscard]] Tensor::Rank1 grad_h(int idx) const;

    // The (symmetric) Cartesian Hessian of h at flat index idx, from the
    // last set_h_derivatives().
    [[nodiscard]] Tensor::Rank2 hess_h(int idx) const;

    // Smallest coordinate distance between neighbouring grid points
    [[nodiscard]] double min_ring_spacing(const std::vector<double> &h) const;

    // Geometric diagnostics

    // Point AHGeometry at AHFinder's persistent surface radius h, physical
    // 3-metric, extrinsic curvature, and its trace. Called once, from
    // AHFinder::init(); all are read through these pointers, so no
    // repeated calls are needed to keep pace with AHFinder's per-step
    // state.
    void set_surface_data(const std::vector<double> *h,
                          const std::vector<Tensor::Rank2> *gamma_LL,
                          const std::vector<Tensor::Rank2> *K_LL,
                          const std::vector<double> *K);

    // Outward unit conormal s_i = F_i / lambda at flat index idx, where
    // F_i = n_i - (grad h)_i is the (un-normalised) level-set gradient of
    // the surface and lambda = sqrt(gamma^ij F_i F_j) normalises it w.r.t.
    // the physical metric. Same construction AHFinder::compute_theta()
    // uses internally for Theta, rebuilt here from direction()/grad_h()
    // rather than shared with it.
    [[nodiscard]] Tensor::Rank1 s_L(int idx) const;

    // Outward unit normal s^i = gamma^ij s_L(idx)_j at flat index idx --
    // the form needed to contract against a lower-index tensor like
    // K_ij or pi_ij = K_ij - K gamma_ij.
    [[nodiscard]] Tensor::Rank1 s_U(int idx) const;

    // pi_ij = K_ij - K gamma_ij at flat index idx: the momentum-density
    // tensor used to build both linear and angular momentum surface
    // integrals -- contracting it with a translational xi^i gives linear
    // momentum, with a rotational xi^i gives angular momentum (spin).
    [[nodiscard]] Tensor::Rank2 pi_LL(int idx) const;

    // Induced 2-metric q_ab of the surface at ring-grid point (i, j), with
    // a, b indexing {theta, phi}: q_ab = gamma_ij e_a^i e_b^j.
    [[nodiscard]] QAB q_ab(int i, int j) const;

    // sqrt(det q_ab) at ring-grid point (i, j).
    [[nodiscard]] amrex::Real root_det_q(int i, int j) const;

    // Total proper area of the surface: the sum of
    // root_det_q(i, j) * dtheta * dphi over the whole ring grid.
    // Not const: refreshes m_dh_dtheta / m_dh_dphi from the currently
    // bound h first, so they can never be read stale.
    amrex::Real area();

    // ADM-type linear momentum of the surface: for each Cartesian
    // direction k, P_k = (1/8pi) * integral pi_kj s^j dA, dA = root_q *
    // dtheta * dphi, using the constant translation vector
    // xi_(k)^i = delta_k^i (so pi_ij xi_(k)^i s^j reduces to pi_kj s^j).
    // Not const, for the same reason as area().
    [[nodiscard]] Tensor::Rank1 momentum();
};

#include "AHGeometry.impl.hpp"

#endif /* AHGEOMETRY_HPP_ */
