/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(PARTICLEINTERPOLATOR_HPP_)
#error "This file should only be included through ParticleInterpolator.hpp"
#endif

#ifndef PARTICLEINTERPOLATOR_IMPL_HPP_
#define PARTICLEINTERPOLATOR_IMPL_HPP_

#include "InterpolationQueryParticle.hpp"
#include "Lagrange.hpp"
#include "StateTypes.hpp"
#include "StateVariables.hpp"
#include "VariableType.hpp"
// amrex includes

#include <AMReX_AmrLevel.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>

// initialise everything and perform some sanity checks
template <int num_reals, int num_components>
void ParticleInterpolator<num_reals, num_components>::setup(GRAmr *gramr_ptr)
{
    // is GRAmr properly set?
    AMREX_ASSERT(gramr_ptr != nullptr);
    m_gr_amr_ptr = gramr_ptr;
    m_bc_params.fill_params();
    GRParmParse particle_interpolator_pp("particle_interpolator");
    particle_interpolator_pp.get("verbosity", m_verbosity);

    this->Define(dynamic_cast<amrex::ParGDBBase *>(m_gr_amr_ptr->GetParGDB()));
    m_initialized = true;

    // read in the physical bounds for reflective BC checks (it is sufficient to
    // do this on lev = 0)
    const amrex::Geometry &geom0 = m_gr_amr_ptr->getLevel(0).Geom();
    m_prob_lo = geom0.RoundOffLo(); // use rounded-off low boundary
    m_prob_hi = geom0.RoundOffHi(); // use rounded-off high boundary

    // set the reflective flags from BC params
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
    {
        m_lo_boundary_reflective[dir] = (m_bc_params.lo_condition[dir] ==
                                         BoundaryConditions::REFLECTIVE_BC);
        m_hi_boundary_reflective[dir] = (m_bc_params.hi_condition[dir] ==
                                         BoundaryConditions::REFLECTIVE_BC);
    }

    // get resolution on level 0
    m_coarsest_dx = geom0.CellSizeArray();
}

// a parity helper (the same way as it was defined in the AMRInterpolator)
template <int num_reals, int num_components>
int ParticleInterpolator<num_reals, num_components>::get_var_parity(
    int comp, int point_idx, const InterpolationQueryParticle &query,
    const Derivative &deriv, VariableType variable_type,
    BCParity derived_parity) const
{
    int parity = 1;

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
    {
        // get the coords
        const amrex::ParticleReal x = query.m_coords[dir][point_idx];

        // check where we are w.r.t to the problem domain
        const bool beyond_lo =
            (m_lo_boundary_reflective[dir] && x < m_prob_lo[dir]);
        const bool beyond_hi =
            (m_hi_boundary_reflective[dir] && x > m_prob_hi[dir]);

        if (beyond_lo || beyond_hi)
        {
            if (variable_type == VariableType::state)
            {
                parity *= BoundaryConditions::get_state_var_parity(comp, dir);
            }
            else if (variable_type == VariableType::derived)
            {
                // if parity was not set for derived vars, print a message
                AMREX_ALWAYS_ASSERT(derived_parity != BCParity::undefined);
                auto dir_parities  = bc_parity_map.at(derived_parity);
                parity            *= dir_parities[dir];
            }
            // invert parity for first derivatives
            if (deriv[dir] == 1)
            {
                parity *= -1;
            }
        }
    }
    return parity;
}

// a function to reflect a particle back into the valid domain, when symmetry
// BCs are used
template <int num_reals, int num_components>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
ParticleInterpolator<num_reals, num_components>::reflect_particle(amrex::Real x,
                                                       amrex::Real low,
                                                       amrex::Real high,
                                                       bool low_reflect,
                                                       bool high_reflect)
{
    // enforce a new particle position if needed
    amrex::Real xl = x;
    if (low_reflect && xl < low)
    {
        xl = low + (low - xl); // reflect across low
    }
    if (high_reflect && xl > high)
    {
        xl = high - (xl - high); // reflect across high
    }

    return xl;
}

// allocate particles at the query points
template <int num_reals, int num_components>
void ParticleInterpolator<num_reals, num_components>::populate_from_query(
    const InterpolationQueryParticle &query)
{
    AMREX_ASSERT(m_initialized);

    m_num_query_points = query.numPoints(); // register number of query points

    // Some ranks can have zero points queried, so return here
    if (query.m_num_points == 0)
    {
        m_particles_populated = true;
        return;
    }

    const int myproc = amrex::ParallelDescriptor::MyProc();

    const int lev  = 0;
    const int grid = 0;
    const int tile = 0;

    // it does not matter on which level the particles are initialised so
    // long as we Redistribute() them later
    auto &ptile = this->DefineAndReturnParticleTile(lev, grid, tile);
    ptile.resize(query.m_num_points);
    auto particle_data = ptile.getParticleTileData();

    // get coords from query
    amrex::GpuArray<const amrex::ParticleReal *, AMREX_SPACEDIM> query_coords{
        query.m_coords[0], query.m_coords[1]
#if (AMREX_SPACEDIM == 3)
        ,
        query.m_coords[2]
#endif
    };

    // Run a check on coords you are interpolating on
    for (int i = 0; i < int(query.m_num_points); ++i)
    {
        amrex::GpuArray<amrex::ParticleReal, AMREX_SPACEDIM> coords;
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
        {
            coords[d] = query_coords[d][i];
        }
        check_domain(coords, 0);
    }

    // copy stuff as otherwise I get errors with lambdas when I use class
    // members inside the GPU loops; this is a bit annoying as I have to
    // copy a few things; any other alternative ways?
    const auto prob_lo    = m_prob_lo;
    const auto prob_hi    = m_prob_hi;
    const auto lo_reflect = m_lo_boundary_reflective;
    const auto hi_reflect = m_hi_boundary_reflective;

    // coords on device
    amrex::GpuArray<amrex::Gpu::DeviceVector<amrex::ParticleReal>,
                    AMREX_SPACEDIM>
        coords_d;
    amrex::GpuArray<const amrex::ParticleReal *, AMREX_SPACEDIM> coords_d_ptr{};

    // copy coords to device vectors
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
    {
        coords_d[d].resize(query.m_num_points);

        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, query_coords[d],
                              query_coords[d] + query.m_num_points,
                              coords_d[d].begin());

        coords_d_ptr[d] = coords_d[d].data(); // get associated device pointer
    }
    amrex::Gpu::streamSynchronize(); // ensure copies complete

    // loop over particles and place them at the required points
    amrex::ParallelFor(
        query.m_num_points,
        [=] AMREX_GPU_DEVICE(int ip)
        {
            auto &p = particle_data[ip];
            p.id()  = ip;     // this is a local index!
            p.cpu() = myproc; // doesn't change on redistribute

            // reflect into valid region and set
            for (int d = 0; d < AMREX_SPACEDIM; ++d)
            {
                p.pos(d) = reflect_particle(
                    static_cast<amrex::Real>(coords_d_ptr[d][ip]), prob_lo[d],
                    prob_hi[d], lo_reflect[d], hi_reflect[d]);
            }

            for (int s = 0; s < num_components; ++s)
            {
                particle_data.rdata(s)[ip] =
                    0.0; // for now set all to zero (this is where we will
                         // store the interpolated values)
            }
        });

    amrex::Gpu::streamSynchronize();

    m_particles_populated = true;
}

// a helper function that helps with interpolation from grid onto particles
template <int num_reals, int num_components>
void ParticleInterpolator<num_reals, num_components>::interpolate_to_particle(
    int lev, amrex::MultiFab &mfab, const amrex::Geometry &geom,
    const InterpolationQueryParticle &query)
{
    const int ncomp = num_components;

    AMREX_ASSERT(num_components == query.numComps());

    if (this->NumberOfParticlesAtLevel(lev) == 0)
        return;

    const auto problem_domain_lo = geom.ProbLoArray();
    const auto dxi               = geom.InvCellSizeArray();
    const auto lo_reflective     = m_lo_boundary_reflective;
    const auto hi_reflective     = m_hi_boundary_reflective;

    // number of cells on the current level
    // (this is needed for handling the higher boundary with symmetric BC in
    // Lagrange interpolation)
    amrex::GpuArray<int, AMREX_SPACEDIM> domain_ncell{};

    for (int d = 0; d < AMREX_SPACEDIM; ++d)
    {
        domain_ncell[d] = geom.Domain().length(d);
    }

    // Gather comp map into managed arrays
    const int num_derivs = static_cast<int>(
        std::distance(query.compsBegin(), query.compsEnd()));

    amrex::Gpu::ManagedVector<Derivative> derivs(num_derivs);
    amrex::Gpu::ManagedVector<int> comp_counts(num_derivs);

    // Flatten comps and generate pointers into it
    amrex::Gpu::ManagedVector<InterpolationQueryParticle::out_t> comps_flat(
        ncomp);
    amrex::Gpu::ManagedVector<InterpolationQueryParticle::out_t *> comps_ptr(
        num_derivs);

    int i = 0, off = 0;
    for (auto comps_it = query.compsBegin(); comps_it != query.compsEnd();
         ++comps_it, ++i)
    {
        derivs[i] = comps_it->first;
        // comp_counts tracks the number of components interpolated
        // for each derivative; e.g. if we have:
        // LOCAL: [chi, K, A_11]
        // DX: [chi]
        // Then comp_counts[0] = 3 and comp_counts[1] = 1
        comp_counts[i] = static_cast<int>(comps_it->second.size());
        comps_ptr[i]   = comps_flat.data() + off;
        for (const auto &entry : comps_it->second)
        {
            comps_flat[off++] = entry;
        }
    }

    Derivative *derivs_data                              = derivs.data();
    InterpolationQueryParticle::out_t *const *comps_data = comps_ptr.data();
    const int *comp_counts_data                          = comp_counts.data();

    amrex::GpuArray<bool, AMREX_SPACEDIM> need_d1{};
    amrex::GpuArray<bool, AMREX_SPACEDIM> need_d2{};
    for (int di = 0; di < num_derivs; ++di)
    {
        for (int dim = 0; dim < AMREX_SPACEDIM; ++dim)
        {
            if (derivs[di][dim] == 1)
            {
                need_d1[dim] = true;
            }
            else if (derivs[di][dim] == 2)
            {
                need_d2[dim] = true;
            }
        }
    }

    // loop over tiles and interpolate now
    for (ParIterType par_iter(*this, lev); par_iter.isValid(); ++par_iter)
    {
        auto &particle_tile     = this->ParticlesAt(lev, par_iter);
        auto particle_tile_data = particle_tile.getParticleTileData();
        const int num_particles = par_iter.numParticles();
        auto fab_array          = mfab[par_iter].const_array();

        amrex::ParallelFor(
            num_particles,
            [=] AMREX_GPU_DEVICE(int ip)
            {
                auto &particle = particle_tile_data[ip];

                amrex::IntVect is_nodal = amrex::IntVect::TheZeroVector();
                // 4th-order Lagrange (5-point stencil)
                Lagrange<s_interp_order + 1>
                    lagrange_interp; // 4th order interpolation
                lagrange_interp.compute_weights(
                    particle, problem_domain_lo, dxi, is_nodal, lo_reflective,
                    hi_reflective, domain_ncell, need_d1, need_d2);

                amrex::ParticleReal interpolated_vals[ncomp];
                lagrange_interp.interpolate(&fab_array, interpolated_vals,
                                            derivs_data, comps_data,
                                            comp_counts_data, num_derivs, dxi);

                // write results to SOA
                for (int icomp = 0; icomp < ncomp; ++icomp)
                {
                    particle_tile_data.rdata(icomp)[ip] =
                        interpolated_vals[icomp];
                }
            });
    }
    // synchronize GPU streams to ensure all particles are updated
    amrex::Gpu::streamSynchronize();
}

// A wrapper function that the user needs to call to perform interpolation
// It uses/collates together all the methods defined in this class
template <int num_reals, int num_components>
void ParticleInterpolator<num_reals, num_components>::interp(
    const InterpolationQueryParticle &query, bool a_refresh_particles,
    const std::string &name_derived, amrex::Real time_derived /*=0.0*/)
{
    AMREX_ASSERT(m_initialized);

    if (query.numComps() > num_components)
    {
        std::string msg =
            "ParticleInterpolator::interp() Oi oi oi! Your query asks for " +
            std::to_string(query.numComps()) +
            " components but this ParticleInterpolator was declared with only "
            "num_components = " +
            std::to_string(num_components) +
            ". Each requested derivative counts as a separate "
            "component.";

        amrex::Abort(msg);
    }

    // Populate particles
    // here we need to cover various scenarious:
    // (1) first call, refresh=false -> populate once + redistribute (automatic
    // since we start with m_need_redistribute=true) (2) first call,
    // refresh=true -> populate once + redistribute (3) later call, refresh=true
    // -> clear + repopulate + redistribute (need to force
    // m_need_redistribute=true)
    if (!m_particles_populated || a_refresh_particles)
    {
        if (a_refresh_particles && m_particles_populated)
        {
            if (m_verbosity)
            {
                amrex::Print() << "ParticleInterpolator: refreshing particles "
                                  "from query\n";
            }

            this->clearParticles();
            force_redistribute(true); // force a redistribution since we have
                                      // cleared the particles
        }
        else if (m_verbosity)
        {
            amrex::Print()
                << "ParticleInterpolator: populating particles from query\n";
        }

        populate_from_query(query);
    }

    AMREX_ALWAYS_ASSERT(query.numPoints() ==
                        m_num_query_points); // check that the number of query
                                             // points has not changed
    ensure_redistributed();

    VariableType variable_type = query.getVariableType();

    // Interpolate to all particles
    if (variable_type == VariableType::state)
    {
        const auto &unique_comps = query.uniqueComps();

        for (int lev = 0; lev <= m_gr_amr_ptr->finestLevel(); ++lev)
        {
            if (this->NumberOfParticlesAtLevel(lev) == 0)
                continue;

            amrex::AmrLevel &level = m_gr_amr_ptr->getLevel(lev);
            amrex::Real cur_time = level.get_state_data(state_index).curTime();
            const amrex::Geometry &geom = level.Geom();
            amrex::MultiFab &state      = level.get_new_data(state_index);

            // Fill ghost cells
            // So FillPatch and FillBoundary in amrex are different routines!
            // FillPatch is more general and fills the ghost cells at the
            // boundaries of fine and coarse levels, whilst FillBoundary is
            // single-level operation only! There is a nice explanation on this
            // issue here: https://github.com/AMReX-Codes/amrex/issues/391
            // We fill one component at a time so that we do not rely on the
            // queried comps being contiguous
            for (const int comp : unique_comps)
            {
                AMREX_ALWAYS_ASSERT(comp >= 0 && comp < state.nComp());
                amrex::AmrLevel::FillPatch(level, state, s_num_ghosts, cur_time,
                                           state_index, comp, 1, comp);
            }

            interpolate_to_particle(lev, state, geom, query);
        }
    }
    // Interpolation for derived vars, takes in MultiFab and comps (unique
    // components numbers), as e.g. we may want to interpolate several fields at
    // once
    else if (variable_type == VariableType::derived)
    {
        // If you look into AMReX_AmrLevel.cpp you will find that derive calls
        // FillPatch automatically, so no need to worry about ghost cells for
        // derived vars.
        auto out_derived =
            m_gr_amr_ptr->derive(name_derived, time_derived, s_num_ghosts);
        amrex::Vector<amrex::MultiFab *> derived_mf_vect;
        // convert vector of unique_ptrs to one of raw pointers
        derived_mf_vect = amrex::GetVecOfPtrs(out_derived);

        for (int lev = 0; lev <= m_gr_amr_ptr->finestLevel(); ++lev)
        {
            if (this->NumberOfParticlesAtLevel(lev) == 0)
                continue;

            auto &level      = m_gr_amr_ptr->getLevel(lev);
            const auto &geom = level.Geom();
            auto &mf         = *derived_mf_vect[lev];

            interpolate_to_particle(lev, mf, geom, query);
        }
    }
    else
    {
        amrex::Abort(
            "ParticleInterpolator::interp(): unsupported VariableType");
    }

    // Aggregate results
    aggregate_points(query);
}

// A function that puts the logic of mpi send and receive buffers together and
// prepares the final out arrays from interpolation
template <int num_reals, int num_components>
void ParticleInterpolator<num_reals, num_components>::aggregate_points(
    const InterpolationQueryParticle &query)
{
    AMREX_ASSERT(m_initialized);

    // pack m_answer_idx and m_answer_data
    prepare_send_buffers();
    // initialise m_query_idx and m_query_data
    prepare_receive_buffers();
    // exchange answers
    exchange_answers();
    // apply parity
    apply_parity_and_store_values(query);
}

// Prepare send buffers and package them up: who do I send answers to?
template <int num_reals, int num_components>
void ParticleInterpolator<num_reals, num_components>::prepare_send_buffers()
{
    const int nprocs = amrex::ParallelDescriptor::NProcs();

    m_mpi.clearAnswerCounts();

    // a temporary storage vector for component values, one entry per local
    // particle
    std::array<std::vector<amrex::ParticleReal>, num_components> comp_values;
    // layout data for storing query ranks and indices
    std::vector<int> query_ranks;
    std::vector<int> query_indices;

    // loop over particles: count + cache data
    for (int lev = 0; lev <= this->finestLevel(); ++lev)
    {
        for (ParIterType par_iter(*this, lev); par_iter.isValid(); ++par_iter)
        {
            const auto &particle_tile = this->ParticlesAt(lev, par_iter);
            const auto particle_aos   = particle_tile.GetArrayOfStructs();
            const int num_particles   = par_iter.numParticles();

            // Copy AoS to host; this is needed as we will need some information
            // (e.g. cpu()) that is stored in the AoS!
            amrex::Gpu::HostVector<ParticleType> particle_aos_h(num_particles);
            amrex::Gpu::copyAsync(amrex::Gpu::deviceToHost, particle_aos.data(),
                                  particle_aos.data() + num_particles,
                                  particle_aos_h.begin());

            // Copy SoA real components to host
            std::array<amrex::Gpu::HostVector<amrex::ParticleReal>,
                       num_components>
                host_soa_real;

            for (int k = 0; k < num_components; ++k)
            {
                host_soa_real[k].resize(num_particles);
                auto soa_real_dptr =
                    particle_tile.GetStructOfArrays().GetRealData(k).data();
                amrex::Gpu::copyAsync(amrex::Gpu::deviceToHost, soa_real_dptr,
                                      soa_real_dptr + num_particles,
                                      host_soa_real[k].begin());
            }

            amrex::Gpu::streamSynchronize();

            const int data_offset = static_cast<int>(
                query_ranks.size()); // offset in the arrays accumulated over
                                     // previous levels/tiles

            // Resize query_ranks and query_indices (and append the space for
            // particles in this level/tile)
            query_ranks.resize(data_offset + num_particles);
            query_indices.resize(data_offset + num_particles);

            // each component vector has num_particles
            for (int k = 0; k < num_components; ++k)
            {
                comp_values[k].resize(data_offset + num_particles);
            }

            for (int i = 0; i < num_particles; ++i)
            {
                const auto &p        = particle_aos_h[i];
                const int query_rank = static_cast<int>(
                    p.cpu()); // this is the rank that has a query, so it
                              // should receive its value

                AMREX_ASSERT(query_rank >= 0 && query_rank < nprocs);

                // how many answers each querying rank will receive
                m_mpi.incrementAnswerCount(query_rank);
                // write in
                const int j    = data_offset + i; // shift the index
                query_ranks[j] = query_rank;
                query_indices[j] =
                    p.id(); // need particle identifier to know which particle
                            // the interpolated value corresponds to

                // cache component values
                for (int k = 0; k < num_components; ++k)
                {
                    comp_values[k][j] =
                        static_cast<amrex::ParticleReal>(host_soa_real[k][i]);
                }
            }

            if (m_verbosity)
            {
                amrex::AllPrint()
                    << "ParticleInterpolator: rank "
                    << amrex::ParallelDescriptor::MyProc() << " holds "
                    << num_particles << " particles on level " << lev << "\n";
            }
        }
    }

    // build a global MPI communication layout: once again answers are sending
    // buffers, queries are receiving buffers
    m_mpi.exchangeLayout();

    const int total_send = m_mpi.totalAnswerCount(); // total answers/sends
    AMREX_ASSERT(static_cast<int>(query_ranks.size()) == total_send);
    AMREX_ASSERT(static_cast<int>(query_indices.size()) == total_send);
    for (int k = 0; k < num_components; ++k)
    {
        AMREX_ASSERT(static_cast<int>(comp_values[k].size()) == total_send);
    }

    m_answer_idx.resize(total_send); // resize the index array I will send
    m_answer_data.resize(
        num_components); // resize the component data array I will send

    // but my m_answer data is in fact
    // m_answer_data[num_components][total_send]
    for (int k = 0; k < num_components; ++k)
    {
        m_answer_data[k].assign(total_send, 0.0); // initialize to zero here
    }

    const int mpi_procs = m_mpi.comm_size();
    std::vector<int> rank_counter(mpi_procs,
                                  0); // to keep track of how many answers I
                                      // have packed for destination rank r

    // Pack the cached data into the flat answer arrays according to MPI layout
    for (int iquery = 0; iquery < total_send; ++iquery)
    {
        const int query_rank = query_ranks[iquery];

        // idx = start of send buffer + how many items I have packaged already
        const int idx =
            m_mpi.answerDispl(query_rank) + rank_counter[query_rank]++;

        m_answer_idx[idx] = query_indices[iquery];

        for (int k = 0; k < num_components; ++k)
        {
            m_answer_data[k][idx] = comp_values[k][iquery];
        }
    }
}

// Prepare receive buffers
template <int num_reals, int num_components>
void ParticleInterpolator<num_reals, num_components>::prepare_receive_buffers()
{
    const int total_recv = m_mpi.totalQueryCount();

    // resize recv buffers (m_query_*) to what THIS rank will receive
    m_query_idx.resize(total_recv);
    m_query_data.resize(num_components);

    for (int k = 0; k < num_components; ++k)
    {
        m_query_data[k].assign(total_recv, 0.0);
    }
}

// Exchange answers between ranks
template <int num_reals, int num_components>
void ParticleInterpolator<num_reals, num_components>::exchange_answers()
{
#ifdef AMREX_USE_MPI
    m_mpi.asyncBegin();

    // exchange indices
    m_mpi.asyncExchangeAnswer(m_answer_idx.data(), m_query_idx.data(), MPI_INT);

    // exchange values for each component
    MPI_Datatype mpi_real =
        amrex::ParallelDescriptor::Mpi_typemap<amrex::Real>::type();

    for (int k = 0; k < num_components; ++k)
    {
        m_mpi.asyncExchangeAnswer(m_answer_data[k].data(),
                                  m_query_data[k].data(), mpi_real);
    }

    m_mpi.asyncEnd();
#else
    // serial
    m_query_idx  = m_answer_idx;
    m_query_data = m_answer_data;
#endif
}

// Build values at particle positions and apply parities
template <int num_reals, int num_components>
void ParticleInterpolator<num_reals, num_components>::apply_parity_and_store_values(
    const InterpolationQueryParticle &query)
{

    const int num_points = static_cast<int>(query.numPoints());
    const int total_recv = m_mpi.totalQueryCount();

    AMREX_ALWAYS_ASSERT(total_recv == num_points);

    // Build a mapping between query index ip and position i in
    // m_query_data[?][i]
    std::vector<int> mpi_mapping; // size of num_points, maps ip to recv index
    mpi_mapping.assign(num_points, -1);

    for (int i = 0; i < total_recv; ++i)
    {
        const int index = m_query_idx[i]; // local query index on this rank

        AMREX_ALWAYS_ASSERT(index >= 0);
        AMREX_ALWAYS_ASSERT(
            index <
            num_points); // because mpi_mapping is of size num_points minus 1!

        mpi_mapping[index] = i;
    }

#if AMREX_DEBUG
    for (int ip = 0; ip < num_points; ++ip)
    {
        if (mpi_mapping[ip] < 0)
        {
            amrex::AllPrint() << "ParticleInterpolator: Rank "
                              << amrex::ParallelDescriptor::MyProc()
                              << " missing answer for ip=" << ip << "\n";
            amrex::Abort(
                "ParticleInterpolator: Missing answers. Mapping incomplete.");
        }
    }
#endif

    // Apply parity

    int comp_idx = 0;
    for (auto deriv_it = query.compsBegin(); deriv_it != query.compsEnd();
         ++deriv_it)
    {
        const auto &comps = deriv_it->second;

        const Derivative &dkey = deriv_it->first;

        const VariableType variable_type = query.getVariableType();

        for (auto &entry : comps)
        {
            const int comp           = entry.comp;
            amrex::ParticleReal *out = entry.out_data_ptr;

            for (int ip = 0; ip < num_points; ++ip)
            {
                const int recv_idx = mpi_mapping[ip];

                AMREX_ASSERT(recv_idx >= 0);

                int parity = get_var_parity(comp, ip, query, dkey,
                                            variable_type, entry.parity);

                out[ip] = parity * m_query_data[comp_idx][recv_idx];
            }
            comp_idx++;
        }
    }
}

// A function to check whether the query point is inside the physical domain
template <int num_reals, int num_components>
void ParticleInterpolator<num_reals, num_components>::check_domain(
    amrex::GpuArray<amrex::ParticleReal, AMREX_SPACEDIM> &x,
    int guard_cells) const
{
    AMREX_ASSERT(guard_cells >= 0);

    for (int d = 0; d < AMREX_SPACEDIM; ++d)
    {

        const amrex::ParticleReal lo_g =
            m_prob_lo[d] - guard_cells * m_coarsest_dx[d];
        const amrex::ParticleReal hi_g =
            m_prob_hi[d] + guard_cells * m_coarsest_dx[d];

        // reflect into the required side if reflective
        amrex::ParticleReal xr = x[d];
        if (m_lo_boundary_reflective[d] && xr < m_prob_lo[d])
        {
            xr = 2.0 * m_prob_lo[d] - xr; // reflect across lo
        }
        else if (m_hi_boundary_reflective[d] && xr > m_prob_hi[d])
        {
            xr = 2.0 * m_prob_hi[d] - xr; // reflect across hi
        }

        bool error =
            false; // error flag to track when we are outside the domain

        // outside on a non-reflective side (beyond physical domain)
        if (!m_lo_boundary_reflective[d] && x[d] <= lo_g)
            error = true;
        if (!m_hi_boundary_reflective[d] && x[d] >= hi_g)
            error = true;

        // still outside physical domain, even after applying reflective
        // conditions
        if (xr < lo_g || xr > hi_g)
            error = true;

        if (error)
        {
            std::string msg =
                "ParticleInterpolator::check_domain() Oi oi oi! "
                "You are trying to access the point at x[" +
                std::to_string(d) + "] = " + std::to_string(x[d]) +
                " and it lies either outside of the AMReX's rounded-off domain "
                "or right at the physical boundary.";

            amrex::Abort(msg);
        }
    }
}

// Ensure that particles are redistributed if needed
template <int num_reals, int num_components>
void ParticleInterpolator<num_reals, num_components>::ensure_redistributed()
{
    const int nlev = m_gr_amr_ptr->finestLevel() + 1;

    // m_last_redistribute_step is empty at the beginning, so resize
    // also if we add or drop a level, it will also need appropsiate resizing
    if (m_last_redistribute_step.size() != nlev)
    {
        m_last_redistribute_step.resize(
            nlev, -1); // put -1s to indicate no redistribute has happened yet

        m_need_redistribute = true;
    }

    // Did a regrid occur since the last redistribution?
    for (int lev = 0; lev < nlev; ++lev)
    {
        int last_regrid_step =
            m_gr_amr_ptr->levelSteps(lev) - m_gr_amr_ptr->levelCount(lev);

        if (last_regrid_step > m_last_redistribute_step[lev])
        {
            m_need_redistribute = true;
            break;
        }
    }

    int need = (m_need_redistribute ? 1 : 0);
    amrex::ParallelDescriptor::ReduceIntMax(
        need); // do we want all ranks to redistribute particles, if one rank
               // requires to do so?
    if (need)
    {
        if (m_verbosity)
        {
            amrex::Print()
                << "ParticleInterpolator: Redistributing particles\n";
        }
        this->Redistribute();

        for (int lev = 0; lev < nlev; ++lev)
        {
            m_last_redistribute_step[lev] =
                m_gr_amr_ptr->levelSteps(lev) - m_gr_amr_ptr->levelCount(lev);
        }

        m_need_redistribute = false;
    }
}

// Option to force Redistribute() flag if needed globally
template <int num_reals, int num_components>
void ParticleInterpolator<num_reals, num_components>::force_redistribute(
    bool flag)
{
    m_need_redistribute = flag;
}

#endif /* PARTICLEINTERPOLATOR_IMPL_HPP_ */
