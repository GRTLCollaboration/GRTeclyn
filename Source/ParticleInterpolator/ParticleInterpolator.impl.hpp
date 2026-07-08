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
template <int num_components>
void ParticleInterpolator<num_components>::setup(
    GRAMR *gramr_ptr, const BoundaryConditions::params_t &a_bc_params,
    bool a_verbosity)
{
    // is GRAMR properly set?
    AMREX_ASSERT(gramr_ptr != nullptr);
    m_gramr_ptr = gramr_ptr;
    m_bc_params = a_bc_params;
    m_verbosity = a_verbosity;

    this->Define(dynamic_cast<amrex::ParGDBBase *>(m_gramr_ptr->GetParGDB()));
    m_initialized = true;

    // read in the physical bounds for reflective BC checks (it is sufficient to
    // do this on lev = 0)
    const amrex::Geometry &geom0 = m_gramr_ptr->getLevel(0).Geom();
    m_prob_lo = geom0.RoundOffLo(); // use rounded-off low boundary
    m_prob_hi = geom0.RoundOffHi(); // use rounded-off high boundary

    const int num_levels = m_gramr_ptr->finestLevel() + 1;

    // Now write in the number of cells on each level (this is needed for
    // handling the higher boundary with symmetric BC in Lagrange interpolation)
    m_domain_ncell.resize(num_levels);

    for (int lev = 0; lev < num_levels; ++lev)
    {
        const amrex::Geometry &geom = m_gramr_ptr->getLevel(lev).Geom();

        for (int d = 0; d < AMREX_SPACEDIM; ++d)
        {
            m_domain_ncell[lev][d] = geom.Domain().length(d);
        }
    }

    // set the reflective flags from BC params
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
    {
        m_lo_boundary_reflective[dir] =
            (m_bc_params.lo_boundary[dir] == BoundaryConditions::REFLECTIVE_BC);
        m_hi_boundary_reflective[dir] =
            (m_bc_params.hi_boundary[dir] == BoundaryConditions::REFLECTIVE_BC);
    }

    // get resolution on level 0
    m_coarsest_dx = geom0.CellSizeArray();
}

// a parity helper (the same way as it was defined in the AMRInterpolator)
template <int num_components>
int ParticleInterpolator<num_components>::get_var_parity(
    int comp, int point_idx, const InterpolationQueryParticle &query,
    const Derivative &deriv, VariableType variable_type,
    BCParity derived_parity) const
{
    int parity = 1;

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
    {
        // get the coords
        const double x = query.m_coords[dir][point_idx];

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
template <int num_components>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
ParticleInterpolator<num_components>::reflect_particle(amrex::Real x,
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
template <int num_components>
void ParticleInterpolator<num_components>::populate_from_query()
{
    AMREX_ASSERT(m_initialized);
    AMREX_ALWAYS_ASSERT(m_query);

    auto &query = *m_query;

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
    amrex::GpuArray<const double *, AMREX_SPACEDIM> query_coords{
        query.m_coords[0], query.m_coords[1]
#if (AMREX_SPACEDIM == 3)
        ,
        query.m_coords[2]
#endif
    };

    // Run a check on coords you are interpolating on
    for (int i = 0; i < int(query.m_num_points); ++i)
    {
        amrex::GpuArray<double, AMREX_SPACEDIM> coords;
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
    amrex::GpuArray<amrex::Gpu::DeviceVector<double>, AMREX_SPACEDIM> coords_d;
    amrex::GpuArray<const double *, AMREX_SPACEDIM> coords_d_ptr{};

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
template <int num_components>
void ParticleInterpolator<num_components>::interpolate_to_particle(
    int lev, amrex::MultiFab &mfab, const amrex::Geometry &geom)
{
    int start_comp  = get_start_comp();
    const int ncomp = num_components;

    AMREX_ASSERT(mfab.nComp() >= start_comp + ncomp);

    if (this->NumberOfParticlesAtLevel(lev) == 0)
        return;

    const auto problem_domain_lo = geom.ProbLoArray();
    const auto dxi               = geom.InvCellSizeArray();
    const auto lo_reflective     = m_lo_boundary_reflective;
    const auto hi_reflective     = m_hi_boundary_reflective;
    const auto domain_ncell      = m_domain_ncell[lev];

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
                lagrange_interp.compute_weights(particle, problem_domain_lo,
                                                dxi, is_nodal, lo_reflective,
                                                hi_reflective, domain_ncell);

                amrex::ParticleReal interpolated_vals[ncomp];
                lagrange_interp.interpolate(&fab_array, interpolated_vals,
                                            start_comp, ncomp);

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
template <int num_components>
void ParticleInterpolator<num_components>::interp(
    InterpolationQueryParticle &query, const std::string &name_derived,
    double time_derived /*=0.0*/)
{
    // Populate particles
    if (!m_particles_populated)
    {
        if (m_verbosity)
        {
            amrex::AllPrint()
                << "ParticleInterpolator: populating particles from query\n";
        }

        // pass the query over
        m_query = &query;

        populate_from_query();
    }

    AMREX_ASSERT(m_initialized);
    ensure_redistributed();

    VariableType variable_type = query.getVariableType();

    // Interpolate to all particles
    if (variable_type == VariableType::state)
    {
        int start_comp  = get_start_comp();
        const int ncomp = num_components;

        for (int lev = 0; lev <= m_gramr_ptr->finestLevel(); ++lev)
        {
            if (this->NumberOfParticlesAtLevel(lev) == 0)
                continue;

            amrex::AmrLevel &level = m_gramr_ptr->getLevel(lev);
            amrex::Real cur_time = level.get_state_data(state_index).curTime();
            const amrex::Geometry &geom = level.Geom();
            amrex::MultiFab &state      = level.get_new_data(state_index);

            // Fill ghost cells
            // So FillPatch and FillBoundary in amrex are different routines!
            // FillPatch is more general and fills the ghost cells at the
            // boundaries of fine and coarse levels, whilst FillBoundary is
            // single-level operation only! There is a nice explanation on this
            // issue here: https://github.com/AMReX-Codes/amrex/issues/391
            amrex::AmrLevel::FillPatch(level, state, s_num_ghosts, cur_time,
                                       state_index, start_comp, ncomp);

            interpolate_to_particle(lev, state, geom);
        }
    }
    // Interpolation for derived vars, takes in MultiFab and comps (unique
    // components numbers), as e.g. we may want to interpolate several fields at
    // once. However, as per current implementation, we can have only contiguous
    // comps
    else if (variable_type == VariableType::derived)
    {
        // If you look into AMReX_AmrLevel.cpp you will find that derive calls
        // FillPatch automatically, so no need to worry about ghost cells for
        // derived vars.
        auto out_derived =
            m_gramr_ptr->derive(name_derived, time_derived, s_num_ghosts);
        amrex::Vector<amrex::MultiFab *> derived_mf_vect;
        // convert vector of unique_ptrs to one of raw pointers
        derived_mf_vect = amrex::GetVecOfPtrs(out_derived);

        for (int lev = 0; lev <= m_gramr_ptr->finestLevel(); ++lev)
        {
            if (this->NumberOfParticlesAtLevel(lev) == 0)
                continue;

            auto &level      = m_gramr_ptr->getLevel(lev);
            const auto &geom = level.Geom();
            auto &mf         = *derived_mf_vect[lev];

            interpolate_to_particle(lev, mf, geom);
        }
    }
    else
    {
        amrex::Abort(
            "ParticleInterpolator::interp(): unsupported VariableType");
    }

    // Aggregate results
    aggregate_points();
}

// A function that puts the logic of mpi send and receive buffers together and
// prepares the final out arrays from interpolation
template <int num_components>
void ParticleInterpolator<num_components>::aggregate_points()
{
    AMREX_ASSERT(m_initialized);
    AMREX_ALWAYS_ASSERT(m_query);

    // pack m_answer_idx and m_answer_data
    prepare_send_buffers();
    // initialise m_query_idx and m_query_data
    prepare_receive_buffers();
    // exchange answers
    exchange_answers();
    // apply parity
    apply_parity_and_store_values();
}

// Prepare send buffers and package them up: who do I send answers to?
template <int num_components>
void ParticleInterpolator<num_components>::prepare_send_buffers()
{
    const int nprocs = amrex::ParallelDescriptor::NProcs();

    m_mpi.clearAnswerCounts();

    // a temporary storage vector for component values, one entry per local
    // particle
    std::array<std::vector<double>, num_components> comp_values;
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

            const int level_size = static_cast<int>(
                query_ranks
                    .size()); // this should count per level, starts at zero and
                              // then adds as we loop through levels

            // Resize query_ranks and query_indices
            query_ranks.resize(level_size + num_particles);
            query_indices.resize(level_size + num_particles);

            // each component vector has num_particles
            for (int k = 0; k < num_components; ++k)
            {
                comp_values[k].resize(level_size + num_particles);
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
                const int j    = level_size + i; // shift the index
                query_ranks[j] = query_rank;
                query_indices[j] =
                    p.id(); // need particle identifier to know which particle
                            // the interpolated value corresponds to

                // cache component values
                for (int k = 0; k < num_components; ++k)
                {
                    comp_values[k][j] =
                        static_cast<double>(host_soa_real[k][i]);
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
template <int num_components>
void ParticleInterpolator<num_components>::prepare_receive_buffers()
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
template <int num_components>
void ParticleInterpolator<num_components>::exchange_answers()
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
template <int num_components>
void ParticleInterpolator<num_components>::apply_parity_and_store_values()
{
    AMREX_ALWAYS_ASSERT(m_query);

    auto &query    = *m_query;
    int start_comp = get_start_comp();

    const int num_points = static_cast<int>(query.numPoints());
    const int total_recv = m_mpi.totalQueryCount();

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
    for (auto deriv_it = query.compsBegin(); deriv_it != query.compsEnd();
         ++deriv_it)
    {
        using comps_t = std::vector<typename InterpolationQueryParticle::out_t>;
        comps_t &comps = deriv_it->second;

        const Derivative &dkey = deriv_it->first;

        const VariableType variable_type = query.getVariableType();

        for (auto &entry : comps)
        {
            const int comp           = entry.comp;
            amrex::ParticleReal *out = entry.out_data_ptr;

            const int k = comp - start_comp; // reindex from 0
            AMREX_ASSERT(k >= 0 && k < num_components);

            for (int ip = 0; ip < num_points; ++ip)
            {
                const int recv_idx = mpi_mapping[ip];

                int parity = get_var_parity(comp, ip, query, dkey,
                                            variable_type, entry.parity);

                out[ip] = parity * m_query_data[k][recv_idx];
            }
        }
    }
}

// A function to check whether the query point is inside the physical domain
template <int num_components>
void ParticleInterpolator<num_components>::check_domain(
    amrex::GpuArray<double, AMREX_SPACEDIM> &x, int guard_cells) const
{
    AMREX_ASSERT(guard_cells >= 0);

    for (int d = 0; d < AMREX_SPACEDIM; ++d)
    {

        const double lo_g = m_prob_lo[d] - guard_cells * m_coarsest_dx[d];
        const double hi_g = m_prob_hi[d] + guard_cells * m_coarsest_dx[d];

        // reflect into the required side if reflective
        double xr = x[d];
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

// A getter function to get the starting component from the query
template <int num_components>
int ParticleInterpolator<num_components>::get_start_comp()
{
    AMREX_ASSERT(m_query);
    auto it         = m_query->compsBegin();
    const auto &vec = it->second;
    int start_comp  = vec.front().comp;

    return start_comp;
}

// Ensure that particles are redistributed if needed
template <int num_components>
void ParticleInterpolator<num_components>::ensure_redistributed()
{
    const int nlev = m_gramr_ptr->finestLevel() + 1;

    // m_last_redistribute_step is empty at the beginning, so resize
    // also if we add or drop a level, it will also need appropsiate resizing
    if (m_last_redistribute_step.size() != nlev)
    {
        m_last_redistribute_step.resize(
            nlev, -1); // put -1s to indicate no redistribute has happened yet
        // upon initialisation this would automatically trigger a regrid
        m_need_redistribute = true;
    }

    // Did a regrid occur since the last redistribution?
    for (int lev = 0; lev < nlev; ++lev)
    {
        int last_regrid_step =
            m_gramr_ptr->levelSteps(lev) - m_gramr_ptr->levelCount(lev);

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
                m_gramr_ptr->levelSteps(lev) - m_gramr_ptr->levelCount(lev);
        }

        m_need_redistribute = false;
    }
}

// Option to force Redistribute() flag if needed globally
template <int num_components>
void ParticleInterpolator<num_components>::force_redistribute(bool flag)
{
    m_need_redistribute = flag;
}

#endif /* PARTICLEINTERPOLATOR_IMPL_HPP_ */
