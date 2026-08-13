#ifndef AHFINDER_HPP_
#define AHFINDER_HPP_

#include <AMReX_Array.H>
#include <AMReX_ParIter.H>
#include <AMReX_Particles.H>
#include <array>

#include "ParticleInterpolator.hpp"

template <int num_components>
class AHFinder : public ParticleInterpolator<num_components>
{
  private:
    int m_num_particles;

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
    // below which the SER ratio is not trusted. Without these the timestep
    // can collapse in one step when theta crosses zero (see update_v()).
    amrex::Real m_dt_shrink;
    amrex::Real m_dt_grow;
    amrex::Real m_theta_floor;

    // Carries the particles' Cartesian coordinates for the one-time initial
    // seeding in populate_from_query(); not used for interpolation itself.
    InterpolationQueryParticle query;

    // All of the per-particle surface state below is indexed by particle id,
    // i.e. by the query index the particle was seeded from. It is deliberately
    // *not* stored in the particle tiles: the surface stencils in
    // ring_gradient() read a particle's four (theta, phi) neighbours, which in
    // general sit in a different box or on a different rank, so no tile-local
    // view of the data is sufficient. The particles themselves therefore carry
    // nothing but their position, which move_radial() pushes back out to them
    // via update_particle_positions().
    std::vector<double> interp_coords_x{};
    std::vector<double> interp_coords_y{};
    std::vector<double> interp_coords_z{};

    // Radius, radial velocity and pseudo-timestep per particle
    std::vector<double> m_h{};
    std::vector<double> m_v{};
    std::vector<double> m_dt{};

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
    InterpolationQueryParticle m_metric_query_state{};
    InterpolationQueryParticle m_metric_query_deriv{};

    std::vector<double> m_theta_vals{};

    // Theta from the previous iteration, for adapting the timestep
    std::vector<double> m_theta_prev{};

    void generate_spherical_query();

    // Returns the indices of the 4 neighbours of particle j on the ring
    // grid, ordered {north, south, east, west}. East/west wrap around
    // longitude; north/south clamp at the poles (returning j itself).
    amrex::GpuArray<int, 4> neighbours(int j) const;
    void move_radial();
    void init_particle_vals();
    void update_v();
    void h_derivs();
    void h_hessian();
    void setup_metric_query();
    void compute_theta();
    double inf_norm(std::vector<double>);

    // Differentiates a scalar field defined on the ring grid (e.g. h, or
    // one of its Cartesian derivatives) with respect to (theta, phi) via a
    // central difference across the given neighbours, then converts to a
    // Cartesian gradient via the spherical Jacobian at (theta, phi, r).
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
        : m_num_particles(num_particles), query(num_particles),
          interp_coords_x(num_particles), interp_coords_y(num_particles),
          interp_coords_z(num_particles), m_h(num_particles),
          m_v(num_particles), m_dt(num_particles), m_center(center),
          m_guess_radius(guess_radius), m_dhdx(num_particles),
          m_dhdy(num_particles), m_dhdz(num_particles), m_d2h_xx(num_particles),
          m_d2h_yy(num_particles), m_d2h_zz(num_particles),
          m_d2h_xy(num_particles), m_d2h_xz(num_particles),
          m_d2h_yz(num_particles), m_metric_query_state(num_particles),
          m_metric_query_deriv(num_particles), m_theta_vals(num_particles),
          m_theta_prev(num_particles)
    {
    }

    void init(GRAMR *gramr_ptr, const BoundaryConditions::params_t &a_bc_params,
              bool a_verbosity);

    void find();
};

#include "AHFinder.impl.hpp"

#endif /* AHFINDER_HPP_ */
