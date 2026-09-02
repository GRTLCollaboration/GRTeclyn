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

#include "AHFinderState.hpp"
#include "ParticleInterpolator.hpp"

template <int num_components>
class AHFinder : public ParticleInterpolator<num_components>
{
  private:
    int m_num_particles;
    int m_n_local;
    int m_start;

    // Ring (latitude x longitude) decomposition of the particles,
    // computed from m_num_particles. m_n_theta * m_n_phi == m_num_particles.
    int m_n_rings   = 0;
    int m_ring_size = 0;

    // Pseudo-timestepping parameters
    amrex::Real m_eta;
    amrex::Real m_c;
    amrex::Real m_tol;
    amrex::Real m_min_dt;
    amrex::Real m_r;

    // Bounds on the per-iteration change in dt, and the magnitude of theta
    // below which the SER ratio is not trusted.
    amrex::Real m_dt_shrink;
    amrex::Real m_dt_grow;
    amrex::Real m_theta_floor;

    // Coords for particleinterpolator query
    std::vector<double> interp_coords_x{};
    std::vector<double> interp_coords_y{};
    std::vector<double> interp_coords_z{};

    // Normalised directions to AHFinder centre for each particle
    std::vector<double> m_dir_x{};
    std::vector<double> m_dir_y{};
    std::vector<double> m_dir_z{};

    // State storing h and v values for all particles. Stored off particles
    // since we need h from other particles to compute its derivative, and
    // this cannot be accessed from another particle if they are not on the
    // same tile
    AHState m_state{};

    // AMReX time integrator for evolution of h and v.
    std::unique_ptr<amrex::TimeIntegrator<AHState>> m_integrator;

    std::array<double, AMREX_SPACEDIM> m_center;
    double m_guess_radius;

    // Gradient and Hessian of h
    std::vector<double> m_dhdx{};
    std::vector<double> m_dhdy{};
    std::vector<double> m_dhdz{};
    std::vector<double> m_d2h_xx{};
    std::vector<double> m_d2h_yy{};
    std::vector<double> m_d2h_zz{};
    std::vector<double> m_d2h_xy{};
    std::vector<double> m_d2h_xz{};
    std::vector<double> m_d2h_yz{};

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

    void generate_spherical_query();

    // Returns the indices of the 4 neighbours of particle j on the ring
    // grid, ordered {north, south, east, west}
    amrex::GpuArray<int, 4> neighbours(int j) const;
    void init_particle_vals();

    // Set particles' coordinates according to their distance from the centre
    void set_particle_positions(const std::vector<double> &h);

    // RHS function to allow amrex time integrator to update h and v
    void compute_rhs(AHState &rhs, AHState &state, amrex::Real time);

    // Update pseudo-timestep based on ratio of improvement of Theta between
    // steps
    amrex::Real update_dt(amrex::Real dt, double theta_old, double theta_new,
                          const std::vector<double> &h) const;

    void h_derivs(const std::vector<double> &h);
    void h_hessian(const std::vector<double> &h);
    void setup_metric_query();
    void compute_theta(const std::vector<double> &h);
    double inf_norm(std::vector<double>);

    // Differentiates a scalar field defined on the ring grid
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE static amrex::GpuArray<amrex::Real, 3>
    ring_gradient(const double *field_ptr, int north_idx, int south_idx,
                  int east_idx, int west_idx, double d_theta, double d_phi,
                  double theta, double phi, double r);

  public:

    using Base         = ParticleInterpolator<num_components>;
    using ParIterType  = typename Base::ParIterType;
    using ParticleType = typename Base::ParticleType;
    using Base::Base;

    AHFinder(int num_particles, std::array<double, AMREX_SPACEDIM> &center,
             double guess_radius = 1.0)
        : m_num_particles(num_particles), m_n_local(local_count(num_particles)),
          m_start(local_start(num_particles)), interp_coords_x(num_particles),
          interp_coords_y(num_particles), interp_coords_z(num_particles),
          m_dir_x(num_particles), m_dir_y(num_particles),
          m_dir_z(num_particles), m_state(std::vector<double>(num_particles),
                                          std::vector<double>(num_particles)),
          m_center(center), m_guess_radius(guess_radius), m_dhdx(num_particles),
          m_dhdy(num_particles), m_dhdz(num_particles), m_d2h_xx(num_particles),
          m_d2h_yy(num_particles), m_d2h_zz(num_particles),
          m_d2h_xy(num_particles), m_d2h_xz(num_particles),
          m_d2h_yz(num_particles), m_metric_query_state(m_n_local),
          m_metric_query_deriv(m_n_local), m_theta_vals(num_particles)
    {
    }

    void init(GRAmr *gramr_ptr);

    void find();
};

#include "AHFinder.impl.hpp"

#endif /* AHFINDER_HPP_ */
