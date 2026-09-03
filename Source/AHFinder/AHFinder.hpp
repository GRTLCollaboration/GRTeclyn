#ifndef AHFINDER_HPP_
#define AHFINDER_HPP_

#include <AMReX_Array.H>
#include <AMReX_ParIter.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Particles.H>
#include <AMReX_TimeIntegrator.H>
#include <algorithm>
#include <array>
#include <memory>

#include "AHFinderParameters.hpp"
#include "AHFinderState.hpp"
#include "AHGeometry.hpp"
#include "ParticleInterpolator.hpp"
#include "Tensor.hpp"

template <int num_components>
class AHFinder : public ParticleInterpolator<num_components>
{
  private:
    int m_num_particles;
    int m_n_local;
    int m_start;

    // Pseudo-timestepping parameters read from the "ah_finder" scope of the
    // input file (eta, c, tolerance, r, cfl_factor). Filled by init().
    ah_finder_params_t m_params{};

    // Smallest permitted pseudo-timestep, bounds on the per-iteration change
    // in dt, and the magnitude of theta below which the SER ratio is not
    // trusted. Not input parameters: these are guard rails on the adaptive
    // timestep rather than knobs to tune per run.
    amrex::Real m_min_dt;
    amrex::Real m_dt_shrink;
    amrex::Real m_dt_grow;
    amrex::Real m_theta_floor;

    // Coords for particleinterpolator query
    std::vector<double> interp_coords_x{};
    std::vector<double> interp_coords_y{};
    std::vector<double> interp_coords_z{};

    // State storing h and v values for all particles. Stored off particles
    // since we need h from other particles to compute its derivative, and
    // this cannot be accessed from another particle if they are not on the
    // same tile
    AHState m_state{};

    // Owns the ring (latitude x longitude) grid: the per-particle
    // directions, the finite-difference stencil, the derivatives of h, and
    // the surface diagnostics (area)
    AHGeometry m_geometry;

    // Physical 3-metric gamma_ij at each particle (flat-indexed as
    // i * m_geometry.ring_size() + j), computed each step in
    // compute_theta(). AHGeometry is given a pointer to this in init(), so
    // it always reads the latest values without a separate copy.
    std::vector<Tensor::Rank2> m_gamma_LL{};

    // AMReX time integrator for evolution of h and v.
    std::unique_ptr<amrex::TimeIntegrator<AHState>> m_integrator;

    // Output arrays for interpolation queries
    std::array<std::vector<double>, 14> m_metric_state{};
    std::array<std::vector<double>, 7> m_metric_dx{};
    std::array<std::vector<double>, 7> m_metric_dy{};
    std::array<std::vector<double>, 7> m_metric_dz{};

    // Split up queries num_components doubles as both the
    // query's flat scratch-array size and the number of contiguous grid
    // comps FillPatch fetches, so a query's total (comp, derivative) entry
    // count can't exceed the simulation's total number of state variables.
    InterpolationQueryParticle m_metric_query_state;
    InterpolationQueryParticle m_metric_query_deriv;

    std::vector<double> m_theta_vals{};

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

    // RHS function to allow amrex time integrator to update h and v
    void compute_rhs(AHState &rhs, AHState &state, amrex::Real time);

    // Update pseudo-timestep based on ratio of improvement of Theta between
    // steps, capped by a CFL condition on the ring-grid spacing
    amrex::Real update_dt(amrex::Real dt, double theta_old, double theta_new,
                          const std::vector<double> &h) const;

    void setup_metric_query();
    void compute_theta(const std::vector<double> &h);
    double inf_norm(std::vector<double>);

  public:

    using Base         = ParticleInterpolator<num_components>;
    using ParIterType  = typename Base::ParIterType;
    using ParticleType = typename Base::ParticleType;
    using Base::Base;

    AHFinder(int num_particles,
             const std::array<double, AMREX_SPACEDIM> &center,
             double guess_radius = 1.0)
        : m_num_particles(num_particles), m_n_local(local_count(num_particles)),
          m_start(local_start(num_particles)), interp_coords_x(num_particles),
          interp_coords_y(num_particles), interp_coords_z(num_particles),
          m_state(std::vector<double>(num_particles),
                  std::vector<double>(num_particles)),
          m_geometry(num_particles, center, guess_radius),
          m_gamma_LL(num_particles), m_metric_query_state(m_n_local),
          m_metric_query_deriv(m_n_local), m_theta_vals(num_particles)
    {
    }

    void init(GRAmr *gramr_ptr);

    void find();
};

#include "AHFinder.impl.hpp"

#endif /* AHFINDER_HPP_ */
