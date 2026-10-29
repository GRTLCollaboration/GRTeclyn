#ifndef AHFINDER_HPP_
#define AHFINDER_HPP_

#include <AMReX_Array.H>
#include <AMReX_ParIter.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Particles.H>
#include <algorithm>
#include <array>
#include <vector>

#include "AHFinderParameters.hpp"
#include "AHGeometry.hpp"
#include "ParticleInterpolator.hpp"
#include "Tensor.hpp"

// LAPACK banded LU solve
extern "C" void dgbsv_(const int *n, const int *kl, const int *ku,
                       const int *nrhs, double *ab, const int *ldab, int *ipiv,
                       double *b, const int *ldb, int *info);

// Number of components (including derivatives) the AHFinder needs to
// interpolate to calculate the expansion
// (values of chi, h_ij, K, A_ij: 14; first derivatives of chi, h_ij: 7 * 3 =
// 21), all fetched in a single query
static constexpr int ah_finder_num_components = 35;

class AHFinder : public ParticleInterpolator<ah_finder_num_components>
{
  private:
    int m_num_particles;
    int m_n_local;
    int m_start;

    // Pseudo-timestepping parameters read from the "ah_finder" scope of the
    // input file (tolerance, r, initial_dt, max_iter, max_stall_iter).
    // Filled by init().
    ah_finder_params_t m_params{};

    // trusted, and the
    // Smallest permitted pseudo-timestep
    amrex::Real m_min_dt{};
    // Largest allowed per-iteration reduction in dt
    amrex::Real m_dt_shrink{};
    // scale sigma of the finite-difference steps used to
    // build the Jacobian.
    amrex::Real m_fd_sigma{};

    // Coords for particleinterpolator query
    std::vector<double> interp_coords_x{};
    std::vector<double> interp_coords_y{};
    std::vector<double> interp_coords_z{};

    // Surface radius h for all particles. Stored off particles since we
    // need h from other particles to compute its derivative, and this cannot
    // be accessed from another particle if they are not on the same tile
    std::vector<double> m_h{};

    // Owns the ring (latitude x longitude) grid: the per-particle
    // directions, the finite-difference stencil, the derivatives of h, and
    // the surface diagnostics (area)
    AHGeometry m_geometry;

    // Physical 3-metric gamma_ij at each particle (flat-indexed as
    // i * m_geometry.ring_size() + j), computed each step in
    // theta_from_metric(). AHGeometry is given a pointer to this in init(), so
    // it always reads the latest values without a separate copy.
    std::vector<Tensor::Rank2> m_gamma_LL{};

    // Output arrays for interpolation queries
    std::array<std::vector<double>, 14> m_metric_state{};
    std::array<std::vector<double>, 7> m_metric_dx{};
    std::array<std::vector<double>, 7> m_metric_dy{};
    std::array<std::vector<double>, 7> m_metric_dz{};

    // Single query for all metric values and derivatives; its total
    // (comp, derivative) entry count must equal ah_finder_num_components.
    InterpolationQueryParticle m_metric_query;

    std::vector<double> m_theta_vals{};

    // Sparsity, colouring and band layout of the Jacobian J = dTheta/dh,

    // Colour of each column (point) j: no row of J contains two columns of
    // the same colour, so all points of one colour can be perturbed at once.
    std::vector<int> m_colour{};
    int m_n_colours{};

    // m_colour_column[i * m_n_colours + c] is the one column of colour c
    // that row i depends on, or -1 if there is none.
    std::vector<int> m_colour_column{};

    // Band ordering: each ring is folded (0, n-1, 1, n-2, ...) so points
    // either side of the periodic phi seam sit close together.
    // m_band_order[k] is the flat point index at band position k and
    // m_band_index is its inverse.
    std::vector<int> m_band_order{};
    std::vector<int> m_band_index{};

    // Lower and upper bandwidths of J in the band ordering, and the leading
    // dimension of its LAPACK band storage (2 * kl + ku + 1).
    int m_kl{};
    int m_ku{};
    int m_ldab{};

    // Whether the state's ghost cells have been filled during the current
    // find(). The grid is frozen while find() runs, so FillPatch is only
    // needed for the first interpolation of each find().
    bool m_fill_ghosts{true};

    static int local_count(int num_particles)
    {
        const int nprocs = amrex::ParallelDescriptor::NProcs();
        const int myproc = amrex::ParallelDescriptor::MyProc();

        return num_particles / nprocs +
               (myproc < num_particles % nprocs ? 1 : 0);
    }

    static int local_start(int num_particles)
    {
        const int nprocs = amrex::ParallelDescriptor::NProcs();
        const int myproc = amrex::ParallelDescriptor::MyProc();

        return myproc * (num_particles / nprocs) +
               std::min(myproc, num_particles % nprocs);
    }

    void init_particle_vals();

    // Set particles' coordinates according to their distance from the centre
    void set_particle_positions(const std::vector<double> &h);

    // Update pseudo-timestep based on ratio of improvement of Theta between
    // steps.
    amrex::Real update_dt(amrex::Real dt, double theta_old,
                          double theta_new) const;

    void setup_metric_query();

    // Interpolate the metric to the current particle positions
    void interp_metric();

    // Theta at every point for surface radius h with no interpolation. Each
    // rank fills its own points and the result is summed across ranks into
    // theta.
    void theta_from_metric(const std::vector<double> &h,
                           std::vector<double> &theta);

    // Build the stencil colouring, colour lookup table and band ordering of
    // the Jacobian, and measure its bandwidths.
    void setup_jacobian();

    // Position of J(row, col) (flat point indices) in the band storage
    [[nodiscard]] int band_offset(int row, int col) const
    {
        const int r = m_band_index[row];
        const int c = m_band_index[col];
        return m_kl + m_ku + r - c + c * m_ldab;
    }

    // Fill jacobian (LAPACK band storage, m_ldab x m_num_particles,
    // column-major) with J = dTheta/dh at m_h by coloured finite
    // differences. Leaves the particles, the interpolated metric and
    // m_gamma_LL at the perturbed positions m_h + eps, so reposition and
    // reinterpolate before evaluating Theta again.
    void assemble_jacobian(std::vector<double> &jacobian);

    // Solve (J + I/dt) delta_h = -Theta for the pseudo-transient
    // continuation step, using m_theta_vals. dgbsv factorises jacobian in
    // place, so it must be reassembled before the next solve.
    bool lapack_solve(std::vector<double> &jacobian, amrex::Real dt,
                      std::vector<double> &delta_h);
    double inf_norm(std::vector<double>);

  public:

    using Base         = ParticleInterpolator<ah_finder_num_components>;
    using ParIterType  = Base::ParIterType;
    using ParticleType = Base::ParticleType;
    using Base::Base;

    AHFinder(int num_particles,
             const std::array<double, AMREX_SPACEDIM> &center,
             double guess_radius = 1.0)
        : m_num_particles(num_particles), m_n_local(local_count(num_particles)),
          m_start(local_start(num_particles)), interp_coords_x(num_particles),
          interp_coords_y(num_particles), interp_coords_z(num_particles),
          m_h(num_particles), m_geometry(num_particles, center, guess_radius),
          m_gamma_LL(num_particles), m_metric_query(m_n_local),
          m_theta_vals(num_particles)
    {
    }

    void init(GRAmr *gramr_ptr);

    void find();
};

#include "AHFinder.impl.hpp"

#endif /* AHFINDER_HPP_ */
