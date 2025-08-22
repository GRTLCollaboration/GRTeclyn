#ifndef PARTICLEINTERPOLATORS_HPP_
#define PARTICLEINTERPOLATORS_HPP_

#include <AMReX_AmrLevel.H>
#include <AMReX_AmrParGDB.H>

#include "BoundaryConditions.hpp"
#include "GRAMR.hpp"
#include "InterpolationQueryParticle.hpp"

// This class interpolates one variable (that may be multi-component) at
// arbitrary coordinates provided via InterpolationQuery, using amrex particles.

class ParticleInterpolators
    : public amrex::ParticleContainer<
          /*NStructReal*/ 0,             // for positions
          /*NStructInt*/ 1,              // particle index
          /*NArrayReal*/ AMREX_SPACEDIM, // SOA slots to store interpolated
                                         // values, cannot have more than
                                         // AMREX_SPACEDIM for one variable
          /*NArrayInt*/ 0>
{
  private:
    GRAMR *m_gr_amr{nullptr};
    bool m_initialized{
        false};          // a guard to make sure we do not uninitialised GRAMR
    int m_start_comp{0}; // first component
    int m_ncomp{1};      // number of components

    // physical domain corners on level 0 for parity logic
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_prob_lo{};
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_prob_hi{};

    // dx on level 0
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_dx{};

    // reflective BC flags per side on the low and high sides
    amrex::GpuArray<bool, AMREX_SPACEDIM> m_lo_boundary_reflective{{false}};
    amrex::GpuArray<bool, AMREX_SPACEDIM> m_hi_boundary_reflective{{false}};

    // copy of BC params
    BoundaryConditions::params_t m_bc_params{};

  public:
    using amrex::ParticleContainer<0, 1, AMREX_SPACEDIM, 0>::ParticleContainer;

    ParticleInterpolators(const BoundaryConditions::params_t &a_bc_params,
                          int a_start_comp, int a_ncomp);

    // initialise everything and perform some sanity checks
    void set_gramr_ptr(GRAMR *gr_amr_ptr);

    // a parity helper (the same way as it was defined in the AMRInterpolator)
    int get_state_var_parity(int comp, int point_idx,
                             const InterpolationQueryParticle &query,
                             const Derivative &deriv) const;

    // a function to reflect a particle back into the valid domain, when
    // symmetry BCs are used
    amrex::Real reflect_particle(amrex::Real x, amrex::Real lo, amrex::Real hi,
                                 bool lo_reflect, bool hi_reflect) const;

    // allocate particles at the query points
    void populate_from_query(const InterpolationQueryParticle &query);

    // interpolate variables into SOA slots
    void interpolate_to_particle();

    // mirror of AMRInterpolator::interp(); assembles all particle data and
    // writes parity * value into the query out arrays
    void interp(InterpolationQueryParticle &query);

    // A function to check whether the query point is inside the physical domain
    template <int dim>
    void check_domain(const std::array<double, dim> &x, int guard_cells = 0);

    // TODO: I have not tested the below yet!!

    // Should I have a function to clear particles at a specific level?
    void clear_level(int lev);

    // Clear particles on all levels
    void clear_all();
};

#include "ParticleInterpolators.impl.hpp"

#endif // PARTICLEINTERPOLATORS_HPP_
