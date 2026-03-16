#ifndef AHFINDER_HPP_
#define AHFINDER_HPP_

#include <AMReX_Array.H>
#include <AMReX_ParIter.H>
#include <AMReX_Particles.H>

#include "ParticleInterpolator.hpp"

template <int num_components>
class AHFinder : public ParticleInterpolator<2, num_components>
{
  private:
    int m_num_particles;

    // Pseudo-timestepping parameters
    amrex::Real m_eta;
    amrex::Real m_c;
    amrex::Real m_tol;
    amrex::Real m_dt;

    InterpolationQueryParticle query;

    std::vector<double> interp_coords_x{};
    std::vector<double> interp_coords_y{};
    std::vector<double> interp_coords_z{};

    std::vector<double> interp_vals{};

    std::array<double, AMREX_SPACEDIM> m_center;

    void generate_spherical_query();
    void move_radial();
    void init_h_v();
    void update_v();
    double inf_norm(std::vector<double>);

  public:

    using Base         = ParticleInterpolator<2, num_components>;
    using ParIterType  = typename Base::ParIterType;
    using ParticleType = typename Base::ParticleType;
    using Base::Base;

    AHFinder(int num_particles, std::array<double, AMREX_SPACEDIM> &center)
        : m_num_particles(num_particles), query(num_particles),
          interp_coords_x(num_particles), interp_coords_y(num_particles),
          interp_coords_z(num_particles), interp_vals(num_particles),
          m_center(center)
    {
    }

    void init(GRAMR *gramr_ptr, const BoundaryConditions::params_t &a_bc_params,
              bool a_verbosity);
};

#include "AHFinder.impl.hpp"

#endif /* AHFINDER_HPP_ */
