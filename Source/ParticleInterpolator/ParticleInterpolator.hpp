/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef PARTICLEINTERPOLATOR_HPP_
#define PARTICLEINTERPOLATOR_HPP_

#include <AMReX_Array.H>
#include <AMReX_ParIter.H>
#include <AMReX_Particles.H>

#include "BCParity.hpp"
#include "BoundaryConditions.hpp"
#include "GRAmr.hpp"
#include "GRAmrLevel.hpp"
#include "InterpolationQueryParticle.hpp"
#include "MPIContextParticle.hpp"

// This class interpolates one variable (that may be multi-component) at
// arbitrary coordinates provided via InterpolationQuery, using amrex particles.
// It can also interpolate first and second derivatives of components.
// Note that any added derivative becomes a new component, so for example
// if you want to interpolate chi, the first derivative of chi and h_11 the
// value of num_components should be 3

template <int num_components>
class ParticleInterpolator
    : public amrex::ParticleContainer<
          /*NStructReal*/ 0,
          /*NStructInt*/ 0,
          /*NArrayReal*/ num_components, // number of SOA slots to store
                                         // interpolated values (assumes
                                         // contiguous storage)
          /*NArrayInt*/ 0>
{

    static_assert(num_components >= 1);

  protected:
    GRAmr *m_gr_amr_ptr{nullptr};

    // A function to check whether the query point is inside the physical domain
    void check_domain(amrex::GpuArray<amrex::ParticleReal, AMREX_SPACEDIM> &x,
                      int guard_cells = 0) const;

  private:
    bool m_initialized{
        false}; // a guard to make sure we do not uninitialised GRAmr

    // physical domain corners on level 0 for parity logic
    amrex::GpuArray<amrex::ParticleReal, AMREX_SPACEDIM> m_prob_lo{};
    amrex::GpuArray<amrex::ParticleReal, AMREX_SPACEDIM> m_prob_hi{};

    // reflective BC flags per side on the low and high sides
    amrex::GpuArray<bool, AMREX_SPACEDIM> m_lo_boundary_reflective{{false}};
    amrex::GpuArray<bool, AMREX_SPACEDIM> m_hi_boundary_reflective{{false}};

    static constexpr int s_interp_order = 4;
    static constexpr int s_num_ghosts   = s_interp_order / 2;

    bool m_verbosity{false};

    bool m_particles_populated{false};
    std::vector<int>
        m_last_redistribute_step{}; // a vector to keep the steps at which
                                    // redistribute happended (this is a vector
                                    // of values stored for all levels)
    bool m_need_redistribute{true};

    // dx on level 0
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_coarsest_dx{};

    // copy of BC params
    BoundaryConditions::params_t m_bc_params;

    std::size_t m_num_query_points{}; // for storing number of query points
                                      // (later used in the minimal check to
                                      // verify that the query has not changed).
    // it is possible that the user uses the same number of points, but the
    // actual query coords are different, but this would require a more complex
    // infrastructre -- to be considered in the future.
    // we also do not really carefully check for ordering of comps and assume
    // the user knows what they are doing -- this will be changed in Adam's PR.

    // mpi stuff
    MPIContextParticle m_mpi;

    std::vector<int> m_answer_idx{}; // indices of the answers (send buffers)
    std::vector<std::vector<amrex::ParticleReal>>
        m_answer_data{}; // send buffers on the answering rank

    std::vector<int> m_query_idx{}; // indices of query (receiving buffers)
    std::vector<std::vector<amrex::ParticleReal>>
        m_query_data{}; // receive buffers on the query rank

    // a parity helper (the same way as it was defined in the AMRInterpolator)
    [[nodiscard]] int
    get_var_parity(int comp, int point_idx,
                   const InterpolationQueryParticle &query,
                   const Derivative &deriv,
                   VariableType variable_type = VariableType::state,
                   BCParity derived_parity    = BCParity::undefined) const;

    // a function to reflect a particle back into the valid domain, when
    // symmetry BCs are used
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static amrex::Real
    reflect_particle(amrex::Real x, amrex::Real low, amrex::Real high,
                     bool low_reflect, bool high_reflect);

    // A helper function that aggregates all the points together from senders
    // and receivers, collects the them into out arrays and applies parity
    void aggregate_points(const InterpolationQueryParticle &query);

    // A helper function to prepare send buffers, packs m_answer_idx and
    // m_answer_data
    void prepare_send_buffers();

    // A helper function to initialise receive buffers, i.e. m_query_idx and
    // m_query_data
    void prepare_receive_buffers();

    // Use m_mpi to exchange m_answer_* and m_query_* objects
    void exchange_answers();

    // Apply parities and store interpolated values in out arrays
    void apply_parity_and_store_values(const InterpolationQueryParticle &query);

  public:

    using Base         = amrex::ParticleContainer<0, 0, num_components, 0>;
    using ParIterType  = typename Base::ParIterType;
    using ParticleType = typename Base::ParticleType;
    using Base::Base;

    ParticleInterpolator() = default; // default constructible

    // initialise everything and perform some sanity checks
    void setup(GRAmr *gramr_ptr);

    // allocate particles at the query points
    void populate_from_query(const InterpolationQueryParticle &query);

    // A helper function that does interpolation from grid onto particles
    void interpolate_to_particle(int lev, amrex::MultiFab &mfab,
                                 const amrex::Geometry &geom,
                                 const InterpolationQueryParticle &query);

    // final interpolation routine exposed to the users
    // a_refresh_particles flag allows to "refresh" the particles, if e.g. the
    // user opts to change the query in the interpolation
    void interp(const InterpolationQueryParticle &query,
                bool a_refresh_particles, const std::string &name_derived = "",
                amrex::Real time_derived = 0.0);

    void ensure_redistributed();

    void force_redistribute(bool flag);
};

#include "ParticleInterpolator.impl.hpp"

#endif // PARTICLEINTERPOLATOR_HPP_
