/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(PARTICLEINTERPOLATOR_HPP_)
#error "This file should only be included through ParticleInterpolator.hpp"
#endif

#ifndef PARTICLEINTERPOLATOR_IMPL_HPP_
#define PARTICLEINTERPOLATOR_IMPL_HPP_

#include "InterpolationLayoutParticle.hpp"
#include "InterpolationQueryParticle.hpp"
#include "Lagrange.hpp"
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
    GRAMR *gr_amr_ptr, const BoundaryConditions::params_t &a_bc_params,
    bool a_verbosity, const DerivedParity *parities)
{
    // is GRAMR properly set?
    AMREX_ASSERT(gr_amr_ptr != nullptr);
    m_gr_amr    = gr_amr_ptr;
    m_bc_params = a_bc_params;
    m_verbosity = a_verbosity;

    this->Define(dynamic_cast<amrex::ParGDBBase *>(m_gr_amr->GetParGDB()));
    m_initialized = true;

    AMREX_ALWAYS_ASSERT(num_components >= 1);

    // read in the physical bounds for reflective BC checks (it is sufficient to
    // do this on lev = 0)
    const amrex::Geometry &geom0 = m_gr_amr->getLevel(0).Geom();
    m_prob_lo                    = geom0.ProbLoArray();
    m_prob_hi                    = geom0.ProbHiArray();

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

    // set derived var parities if provided
    if (parities != nullptr)
    {
        for (int i = 0; i < num_components; ++i)
        {
            set_derived_var_parity(parities[i].comp, parities[i].parity);
        }
    }
}

template <int num_components>
void ParticleInterpolator<num_components>::set_derived_var_parity(int comp,
                                                                  BCParity p)
{
    AMREX_ALWAYS_ASSERT(comp >= 0 && comp < num_components);
    m_derived_bc_parity[comp] = p;
}

// a parity helper (the same way as it was defined in the AMRInterpolator)
template <int num_components>
int ParticleInterpolator<num_components>::get_var_parity(
    int comp, int point_idx, const InterpolationQueryParticle &query,
    const Derivative &deriv, VariableType variable_type) const
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
                AMREX_ALWAYS_ASSERT(m_derived_bc_parity[comp] !=
                                    BCParity::undefined);
                BCParity comp_parity  = m_derived_bc_parity[comp];
                auto dir_parities     = bc_parity_map.at(comp_parity);
                parity               *= dir_parities[dir];
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
                                                       amrex::Real lo,
                                                       amrex::Real hi,
                                                       bool lo_reflect,
                                                       bool hi_reflect)
{
    // enforce a new particle position if needed
    amrex::Real xl = x;
    if (lo_reflect && xl < lo)
    {
        xl = lo + (lo - xl); // reflect across lo
    }
    if (hi_reflect && xl > hi)
    {
        xl = hi - (xl - hi); // reflect across hi
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

    const int myproc = amrex::ParallelDescriptor::MyProc();

    const int lev = 0;
    const int n =
        static_cast<int>(query.m_num_points); // number of query points

    amrex::MFIter mfi = this->MakeMFIter(lev);
    const int grid    = mfi.index();
    const int tile    = mfi.LocalTileIndex();

    // it does not matter on which level the particles are initialised so
    // long as we Redistribute() them later
    auto &ptile = this->DefineAndReturnParticleTile(lev, grid, tile);
    ptile.resize(n);
    auto particle_data = ptile.getParticleTileData();

    // get coords from query
    const double *x = query.m_coords[0];
    const double *y = query.m_coords[1];
#if AMREX_SPACEDIM == 3
    const double *z = query.m_coords[2];
#endif

    // Run a check on coords you are interpolating on
    for (int i = 0; i < n; i++)
    {
        std::array<double, AMREX_SPACEDIM> coords{
            AMREX_D_DECL(x[i], y[i], z[i])};
        check_domain(coords, 0);
    }

    // copy stuff as otherwise I get errors with lambdas when I use class
    // members inside the GPU loops; this is a bit annoying as I have to
    // copy a few things; any other alternative ways?
    const auto prob_lo    = m_prob_lo;
    const auto prob_hi    = m_prob_hi;
    const auto lo_reflect = m_lo_boundary_reflective;
    const auto hi_reflect = m_hi_boundary_reflective;

    // copy coords here
    const amrex::Real *x_p = nullptr, *y_p = nullptr, *z_p = nullptr;
    amrex::Gpu::DeviceVector<amrex::Real> Xd, Yd, Zd;

    // copy x,y,z to device vectors now
    amrex::Gpu::DeviceVector<double> x_d(n), y_d(n);
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, x, x + n, x_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, y, y + n, y_d.begin());

#if AMREX_SPACEDIM == 3
    amrex::Gpu::DeviceVector<double> z_d(n);
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, z, z + n, z_d.begin());
#endif
    amrex::Gpu::streamSynchronize(); // ensure copies complete

    // Get device pointers to capture by value
    auto x_d_ptr = x_d.data();
    auto y_d_ptr = y_d.data();
#if AMREX_SPACEDIM == 3
    auto z_d_ptr = z_d.data();
#endif

    // loop over particles and place them at the required points
    amrex::ParallelFor(
        n,
        [=] AMREX_GPU_DEVICE(int ip)
        {
            auto &p = particle_data[ip];
            p.id() = ip + 1; // this is a local index! particle id starts from 1
            p.cpu() =
                myproc; // rank number, stays unique "at borth" of the particle

            // reflect into valid region and set
            p.pos(0) = reflect_particle(static_cast<amrex::Real>(x_d_ptr[ip]),
                                        prob_lo[0], prob_hi[0], lo_reflect[0],
                                        hi_reflect[0]);
            p.pos(1) = reflect_particle(static_cast<amrex::Real>(y_d_ptr[ip]),
                                        prob_lo[1], prob_hi[1], lo_reflect[1],
                                        hi_reflect[1]);
#if AMREX_SPACEDIM == 3
            p.pos(2) = reflect_particle(static_cast<amrex::Real>(z_d_ptr[ip]),
                                        prob_lo[2], prob_hi[2], lo_reflect[2],
                                        hi_reflect[2]);
#endif
            p.idata(0) = ip;

            for (int s = 0; s < num_components; ++s)
            {
                particle_data.rdata(s)[ip] =
                    0.0; // for now set all to zero (this is where we will
                         // store the interpolated values)
            }
        });

    amrex::Gpu::streamSynchronize();

    m_particles_seeded = true;

    // amrex::ParallelDescriptor::Barrier(); // TODO! test if this makes a
    // difference.
}

// a helper function that helps with interpolation from grid onto particles
template <int num_components>
void ParticleInterpolator<num_components>::interpolation_to_particle_helper(
    int lev, amrex::MultiFab &mf, const amrex::Geometry &geom, int num_ghosts)
{
    int start_comp                    = start_comp_getter();
    const int ncomp                   = num_components;
    static constexpr int interp_order = 4; // 4th order Lagrange

    AMREX_ASSERT(mf.nComp() >= start_comp + ncomp);

    if (this->NumberOfParticlesAtLevel(lev) == 0)
        return;

    // Fill ghost cells
    amrex::IntVect nghost(AMREX_D_DECL(num_ghosts, num_ghosts, num_ghosts));
    mf.FillBoundary(start_comp, ncomp, nghost, geom.periodicity());

    const auto problem_domain_lo = geom.ProbLoArray();
    const auto dxi               = geom.InvCellSizeArray();

    // loop over tiles and interpolate now
    for (ParIterType par_iter(*this, lev); par_iter.isValid(); ++par_iter)
    {
        auto &particle_tile     = this->ParticlesAt(lev, par_iter);
        auto particle_tile_data = particle_tile.getParticleTileData();
        const int num_particles = par_iter.numParticles();
        auto fab_array          = mf[par_iter].const_array();

        amrex::ParallelFor(
            num_particles,
            [=] AMREX_GPU_DEVICE(int ip)
            {
                auto &particle = particle_tile_data[ip];

                amrex::IntVect is_nodal = amrex::IntVect::TheZeroVector();
                // 4th-order Lagrange (5-point stencil)
                Lagrange<interp_order + 1>
                    lagrange_interp; // 4th order interpolation
                lagrange_interp.compute_weights(particle, problem_domain_lo,
                                                dxi, is_nodal);

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

        // synchronize GPU streams to ensure all particles are updated
        amrex::Gpu::streamSynchronize();
    }
}

// interpolate variables into SOA slots
template <int num_components>
void ParticleInterpolator<num_components>::interpolate_to_particle()
{
    AMREX_ASSERT(m_initialized);

    ensure_redistributed();

    static constexpr int interp_order = 4; // 4th order Lagrange
    static constexpr int num_ghosts =
        interp_order / 2; // number of ghosts needed

    for (int lev = 0; lev <= m_gr_amr->finestLevel(); ++lev)
    {
        if (this->NumberOfParticlesAtLevel(lev) == 0)
            continue;

        amrex::AmrLevel &level      = m_gr_amr->getLevel(lev);
        const amrex::Geometry &geom = level.Geom();
        amrex::MultiFab &state      = level.get_new_data(0);

        interpolation_to_particle_helper(lev, state, geom, num_ghosts);
    }

    m_need_redistribute = false;
}

// Interpolation for derived vars, takes in MultiFab and comps (unique
// components numbers), as e.g. we may want to interpolate several fields at
// once. However, as per current implementation, we can have only contiguous
// comps
template <int num_components>
void ParticleInterpolator<num_components>::
    interpolate_to_particle_from_derived_fields(
        const std::vector<amrex::MultiFab *> &a_derived_mf_vect)
{
    AMREX_ASSERT(m_initialized);

    const int nlevs = m_gr_amr->finestLevel() + 1;
    AMREX_ASSERT((int)a_derived_mf_vect.size() == nlevs);

    ensure_redistributed();

    static constexpr int interp_order = 4; // 4th order Lagrange
    static constexpr int num_ghosts =
        interp_order / 2; // number of ghosts needed

    for (int lev = 0; lev <= m_gr_amr->finestLevel(); ++lev)
    {
        if (this->NumberOfParticlesAtLevel(lev) == 0)
            continue;

        auto &level      = m_gr_amr->getLevel(lev);
        const auto &geom = level.Geom();
        auto &mf         = *a_derived_mf_vect[lev];

        interpolation_to_particle_helper(lev, mf, geom, num_ghosts);
    }

    m_need_redistribute = false;
}

// A wrapper function that the user needs to call to perform interpolation
// It uses/collates together all the methods defined in this class
template <int num_components>
void ParticleInterpolator<num_components>::interp(
    InterpolationQueryParticle &query, VariableType variable_type,
    const std::string &name_derived, double time_derived /*=0.0*/)
{
    static constexpr int interp_order = 4; // 4th order Lagrange
    static constexpr int num_ghosts =
        interp_order / 2; // number of ghosts needed

    // Populate particles
    if (!m_particles_seeded)
    {
        if (m_verbosity)
        {
            amrex::AllPrint()
                << "ParticleInterpolator: populating particles from query\n";
        }
        m_query = std::make_unique<InterpolationQueryParticle>(query);

        populate_from_query();
    }

    // Interpolate to all particles
    if (variable_type == VariableType::state)
    {
        interpolate_to_particle();
    }
    else if (VariableType::derived == variable_type)
    {
        auto out_derived =
            m_gr_amr->derive(name_derived, time_derived, num_ghosts);
        amrex::Vector<amrex::MultiFab *> derived_mf_vect;
        derived_mf_vect = amrex::GetVecOfPtrs(
            out_derived); // see here for more details
                          // (https://amrex-codes.github.io/amrex/doxygen/AMReX__Vector_8H_source.html)

        interpolate_to_particle_from_derived_fields(derived_mf_vect);
    }
    else
    {
        amrex::Abort("The type of VaribaleType is not recognized in "
                     "ParticleInterpolator::interp()!");
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
    send_buffers();
    // initialise m_query_idx and m_query_data
    prepare_receive_buffers();
    // exchange answers
    exchange_answers();
    // build query values and apply parity
    build_values_and_apply_parity();
}

// Prepare send buffers and package them up: who do I send answers to?
template <int num_components>
void ParticleInterpolator<num_components>::send_buffers()
{
    const int nprocs = amrex::ParallelDescriptor::NProcs();

    m_mpi.clearQueryCounts();
    int local_particle_counter = 0;

    // we do not know how many local particles we will see, so start empty.
    InterpolationLayoutParticle layout(0);

    // a temporary storage vector for component values, one entry per local
    // particle
    std::array<std::vector<double>, num_components> comp_values;

    // loop over particles: count + cache data
    for (int lev = 0; lev <= this->finestLevel(); ++lev)
    {
        for (ParIterType par_iter(*this, lev); par_iter.isValid(); ++par_iter)
        {
            const auto &ptile = this->ParticlesAt(lev, par_iter);
            const auto aos    = ptile.GetArrayOfStructs();
            const int np      = par_iter.numParticles();

            if (np == 0)
            {
                continue;
            }

            // Copy AoS to host; this is needed as we will need some information
            // (e.g. cpu()) that is stored in the AoS!
            amrex::Gpu::HostVector<ParticleType> host_aos(np);
            amrex::Gpu::copyAsync(amrex::Gpu::deviceToHost, aos.data(),
                                  aos.data() + np, host_aos.begin());

            // Copy SoA real components to host
            std::array<amrex::Gpu::HostVector<amrex::ParticleReal>,
                       num_components>
                host_soa_real;

            for (int k = 0; k < num_components; ++k)
            {
                host_soa_real[k].resize(np);
                auto soa_real_dptr =
                    ptile.GetStructOfArrays().GetRealData(k).data();
                amrex::Gpu::copyAsync(amrex::Gpu::deviceToHost, soa_real_dptr,
                                      soa_real_dptr + np,
                                      host_soa_real[k].begin());
            }

            amrex::Gpu::streamSynchronize();

            const int myproc = amrex::ParallelDescriptor::MyProc();

            for (int i = 0; i < np; ++i)
            {
                const auto &p        = host_aos[i];
                const int owner_rank = static_cast<int>(
                    p.cpu()); // this is the rank that owns the query, so it
                              // should receive its value
                const int ip = p.idata(0);

                AMREX_ASSERT(owner_rank >= 0 && owner_rank < nprocs);

                // how many answers each owner rank will receive
                m_mpi.incrementQueryCount(owner_rank);
                ++local_particle_counter; // count the particle in
                // cache owner rank and local query index
                layout.rank.push_back(owner_rank);
                layout.q_local.push_back(p.idata(0)); // index in owner's query

                // cache component values
                for (int k = 0; k < num_components; ++k)
                {
                    comp_values[k].push_back(
                        static_cast<double>(host_soa_real[k][i]));
                }
            }
        }
    }

    if (m_verbosity)
    {
        amrex::AllPrint() << "ParticleInterpolator: rank "
                          << amrex::ParallelDescriptor::MyProc() << " holds "
                          << local_particle_counter << " particles\n";
    }

    // build a global MPI communication layout
    m_mpi.exchangeLayout();

    const int total_send =
        m_mpi.totalQueryCount(); // total answers this rank will SEND

    m_answer_idx.resize(total_send); // resize the index array I will send back
    m_answer_data.resize(
        num_components); // resize the component data array I will send back

    // but my m_answer data is in fact
    // m_answer_data[num_components][total_send]
    for (int k = 0; k < num_components; ++k)
    {
        m_answer_data[k].assign(total_send, 0.0); // initialize to zero here
    }

    const int mpi_procs = MPIContext::comm_size();
    std::vector<int> rank_counter(mpi_procs,
                                  0); // to keep track of how many answers I
                                      // have packed for destination rank r

    const int owner_size = static_cast<int>(layout.rank.size());

    // Pack the cached data into the flat answer arrays according to MPI layout
    for (int owner = 0; owner < owner_size; ++owner)
    {
        const int owner_rank = layout.rank[owner];

        // idx = start of send buffer + how many items I have packaged already
        const int idx =
            m_mpi.queryDispl(owner_rank) + rank_counter[owner_rank]++;

        m_answer_idx[idx] = layout.q_local[owner];

        for (int k = 0; k < num_components; ++k)
        {
            m_answer_data[k][idx] = comp_values[k][owner];
        }
    }
}

// Prepare receive buffers
template <int num_components>
void ParticleInterpolator<num_components>::prepare_receive_buffers()
{
    const int total_recv = m_mpi.totalAnswerCount();

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
    m_mpi.asyncExchangeQuery(m_answer_idx.data(), m_query_idx.data(), MPI_INT);

    // exchange values for each component
    MPI_Datatype mpi_real =
        amrex::ParallelDescriptor::Mpi_typemap<amrex::Real>::type();

    for (int k = 0; k < num_components; ++k)
    {
        m_mpi.asyncExchangeQuery(m_answer_data[k].data(),
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
void ParticleInterpolator<num_components>::build_values_and_apply_parity()
{
    AMREX_ALWAYS_ASSERT(m_query);

    auto &query    = *m_query;
    int start_comp = start_comp_getter();

    const int num_points = static_cast<int>(query.numPoints());
    const int total_recv = m_mpi.totalAnswerCount();

    // Build a mapping between query index ip and position i in
    // m_query_data[?][i]
    m_mpi_mapping.assign(num_points, -1);

    for (int i = 0; i < total_recv; ++i)
    {
        const int index = m_query_idx[i]; // local query index on this rank
        if (index < 0 || index >= num_points)
        {
            amrex::Abort(
                "build_values_and_apply_parity(): received index out of range");
        }

        m_mpi_mapping[index] = i;
    }

#if AMREX_DEBUG
    for (int ip = 0; ip < num_points; ++ip)
    {
        if (m_mpi_mapping[ip] < 0)
        {
            amrex::AllPrint() << "Rank " << amrex::ParallelDescriptor::MyProc()
                              << " missing answer for ip=" << ip << "\n";
            amrex::Abort("Missing answers: mapping incomplete");
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

        for (auto &entry : comps)
        {
            const int comp           = entry.comp;
            amrex::ParticleReal *out = entry.out_data_ptr;

            const VariableType variable_type = entry.variable_type;

            const int k = comp - start_comp; // reindex from 0
            AMREX_ASSERT(k >= 0 && k < num_components);

            for (int ip = 0; ip < num_points; ++ip)
            {
                const int recv_idx = m_mpi_mapping[ip];

                int parity =
                    get_var_parity(comp, ip, query, dkey, variable_type);

                out[ip] = parity * m_query_data[k][recv_idx];
            }
        }
    }
}

// A function to check whether the query point is inside the physical domain
template <int num_components>
void ParticleInterpolator<num_components>::check_domain(
    std::array<double, AMREX_SPACEDIM> &x, int guard_cells) const
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
        if (!m_lo_boundary_reflective[d] && x[d] < lo_g)
            error = true;
        if (!m_hi_boundary_reflective[d] && x[d] > hi_g)
            error = true;

        // still outside physical domain, even after applying reflective
        // conditions
        if (xr < lo_g || xr > hi_g)
            error = true;

        if (error)
        {
            std::string msg = "ParticleInterpolator::check_domain() Oi oi oi! "
                              "You are trying to access the point at x[" +
                              std::to_string(d) +
                              "] = " + std::to_string(x[d]) +
                              " and it lies outside of your domain.";

            amrex::Abort(msg);
        }
    }
}

// A getter function to get the starting component from the query
template <int num_components>
int ParticleInterpolator<num_components>::start_comp_getter()
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
    int need = (m_need_redistribute ? 1 : 0);
    amrex::ParallelDescriptor::ReduceIntMax(
        need); // do we want all ranks to redistribute particles, if one rank
               // requires to do so?
    if (need)
    {
        if (m_verbosity)
        {
            amrex::Print() << "Redistributing all particles \n";
        }
        this->Redistribute();
        m_need_redistribute = false;
    }
}

// Option to force Redistribute() flag if needed globally
template <int num_components>
void ParticleInterpolator<num_components>::force_redistribute(bool flag)
{
    m_need_redistribute = flag;
}

// TODO: I have not tested the below yet!!

// Should I have a function to clear particles at a specific level?
template <int num_components>
void ParticleInterpolator<num_components>::clear_level(int lev)
{
    for (ParIterType it(*this, lev); it.isValid(); ++it)
    {
        auto &ptile = this->ParticlesAt(lev, it);
        ptile.resize(0); // by resizing to 0 we clear the particles
    }
    this->Redistribute();
}

// Clear particles on all levels
template <int num_components>
void ParticleInterpolator<num_components>::clear_all()
{
    for (int lev = 0; lev <= this->finestLevel(); ++lev)
    {
        clear_level(lev);
    }
}

#endif /* PARTICLEINTERPOLATOR_IMPL_HPP_ */
