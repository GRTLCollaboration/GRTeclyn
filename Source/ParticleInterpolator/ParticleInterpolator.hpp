#ifndef PARTICLEINTERPOLATOR_HPP_
#define PARTICLEINTERPOLATOR_HPP_

#include <AMReX_Array.H>
#include <AMReX_ParIter.H>
#include <AMReX_Particles.H>

#include "BCParity.hpp"
#include "BoundaryConditions.hpp"
#include "GRAMR.hpp"
#include "InterpolationQueryParticle.hpp"
#include "MPIContext.hpp"

// This class interpolates one variable (that may be multi-component) at
// arbitrary coordinates provided via InterpolationQuery, using amrex particles.

template <int num_components>
class ParticleInterpolator
    : public amrex::ParticleContainer<
          /*NStructReal*/ 0,
          /*NStructInt*/ 1,              // particle index
          /*NArrayReal*/ num_components, // number of SOA slots to store
                                         // interpolated values (assumes
                                         // contiguous storage)
          /*NArrayInt*/ 0>
// TODO: to decide whether we want to store particle index in SoA integer slot
// (see my comments on review page)
{
  private:
    GRAMR *m_gr_amr{nullptr};
    bool m_initialized{
        false}; // a guard to make sure we do not uninitialised GRAMR

    // physical domain corners on level 0 for parity logic
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_prob_lo{};
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_prob_hi{};

    // reflective BC flags per side on the low and high sides
    amrex::GpuArray<bool, AMREX_SPACEDIM> m_lo_boundary_reflective{{false}};
    amrex::GpuArray<bool, AMREX_SPACEDIM> m_hi_boundary_reflective{{false}};

    bool m_verbosity{false};

    bool m_particles_seeded{false};
    bool m_need_redistribute{true};

    std::array<BCParity, num_components> m_derived_bc_parity{
        BCParity::undefined}; // default parity set for derived vars

    // dx on level 0
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_coarsest_dx{};

    // copy of BC params
    BoundaryConditions::params_t m_bc_params{};

    // store the query here
    std::unique_ptr<InterpolationQueryParticle> m_query;
    // for getting the starting component of query
    int start_comp_getter();

    // mpi stuff
    MPIContext m_mpi;
    std::vector<int> m_mpi_mapping; // size of num_points, maps ip to recv index

    std::vector<int> m_answer_idx; // indices of the answers
    std::vector<std::vector<double>>
        m_answer_data; // data buffers on the answering rank

    std::vector<int> m_query_idx; // indices of query
    std::vector<std::vector<double>>
        m_query_data; // receive buffers on the query rank

    // a parity helper (the same way as it was defined in the AMRInterpolator)
    int get_var_parity(int comp, int point_idx,
                       const InterpolationQueryParticle &query,
                       const Derivative &deriv,
                       VariableType variable_type = VariableType::state) const;

    // a function to reflect a particle back into the valid domain, when
    // symmetry BCs are used
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static amrex::Real
    reflect_particle(amrex::Real x, amrex::Real lo, amrex::Real hi,
                     bool lo_reflect, bool hi_reflect);

    // A function to check whether the query point is inside the physical domain
    void check_domain(std::array<double, AMREX_SPACEDIM> &x,
                      int guard_cells = 0) const;

    // helper function to set parities of derived vars per component
    void set_derived_var_parity(int comp, BCParity p);

    // allocate particles at the query points
    void populate_from_query();

    // A helper function that does interpolation from grid onto particles
    void interpolation_to_particle_helper(int lev, amrex::MultiFab &mf,
                                          const amrex::Geometry &geom,
                                          int num_ghosts);

    // Interpolate to particle function
    void interpolate_to_particle();

    // Interpolate to particle for derived variables
    void interpolate_to_particle_from_derived_fields(
        const std::vector<amrex::MultiFab *> &a_derived_mf_vect);

    // A helper function that aggregates all the points together from senders
    // and receivers, collects the them into out arrays and applies parity
    void aggregate_points();

    // A helper function to prepare send buffers, packs m_answer_idx and
    // m_answer_data
    void send_buffers();

    // A helper function to initialise receive buffers, i.e. m_query_idx and
    // m_query_data
    void prepare_receive_buffers();

    // Use m_mpi to exchange m_answer_* and m_query_* objects
    void exchange_answers();

    // Build query values on owner ranks and apply parity into out arrays
    void build_values_and_apply_parity();

  public:

    using Base         = amrex::ParticleContainer<0, 1, num_components, 0>;
    using ParIterType  = typename Base::ParIterType;
    using ParticleType = typename Base::ParticleType;
    using Base::Base;

    ParticleInterpolator() = default; // default constructible

    // A struct to set parity for derived variables
    struct DerivedParity
    {
        int comp;        // component of the derived varable
        BCParity parity; // parity to be applied
    };

    // initialise everything and perform some sanity checks
    void setup(GRAMR *gr_amr_ptr,
               const BoundaryConditions::params_t &a_bc_params,
               bool a_verbosity              = false,
               const DerivedParity *parities = nullptr);

    // final interpolation routine exposed to the users
    void interp(InterpolationQueryParticle &query, VariableType variable_type,
                const std::string &name_derived = "",
                double time_derived             = 0.0);

    void ensure_redistributed();

    void force_redistribute(bool flag);

    // TODO: I have not tested the below yet!!

    // Should I have a function to clear particles at a specific level?
    void clear_level(int lev);

    // Clear particles on all levels
    void clear_all();
};

#include "ParticleInterpolator.impl.hpp"

#endif // PARTICLEINTERPOLATOR_HPP_
