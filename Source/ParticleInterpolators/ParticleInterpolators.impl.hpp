/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(PARTICLEINTERPOLATORS_HPP_)
#error "This file should only be included through ParticleInterpolators.hpp"
#endif

#ifndef PARTICLEINTERPOLATORS_IMPL_HPP_
#define PARTICLEINTERPOLATORS_IMPL_HPP_

#include "InterpolationQueryParticle.hpp"
#include "LagrangeInterpolation.hpp"

// amrex includes

#include <AMReX_AmrLevel.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_ParallelDescriptor.H>

inline ParticleInterpolators::ParticleInterpolators(
    const BoundaryConditions::params_t &a_bc_params, int a_start_comp,
    int a_ncomp)
    : m_gr_amr(nullptr), m_initialized(false), m_start_comp(a_start_comp),
      m_ncomp(a_ncomp), m_bc_params(a_bc_params)
{
    // constructor body
}

// initialise everything and perform some sanity checks
void ParticleInterpolators::set_gramr_ptr(GRAMR *gr_amr_ptr)
{
    // is GRAMR properly set?
    AMREX_ASSERT(gr_amr_ptr != nullptr);
    m_gr_amr = gr_amr_ptr;

    this->Define(dynamic_cast<amrex::ParGDBBase *>(m_gr_amr->GetParGDB()));
    m_initialized = true;

    AMREX_ALWAYS_ASSERT(m_ncomp >= 1 && m_ncomp <= AMREX_SPACEDIM);

    // read in the physical bounds for reflective BC checks (it is sufficient to
    // do this on lev = 0)
    const amrex::Geometry &geom0 = m_gr_amr->getLevel(0).Geom();
    const auto plo               = geom0.ProbLoArray();
    const auto phi               = geom0.ProbHiArray();
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
    {
        m_prob_lo[d] = plo[d];
        m_prob_hi[d] = phi[d];
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
    m_dx = geom0.CellSizeArray();
}

// a parity helper (the same way as it was defined in the AMRInterpolator)
int ParticleInterpolators::get_state_var_parity(
    int comp, int point_idx, const InterpolationQueryParticle &query,
    const Derivative &deriv) const
{
    int parity = 1;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
    {
        // get the coords
        const double x = query.m_coords[dir][point_idx];

        // check where we are w.r.t to the prob domain
        const bool beyond_lo =
            (m_lo_boundary_reflective[dir] && x < m_prob_lo[dir]);
        const bool beyond_hi =
            (m_hi_boundary_reflective[dir] && x > m_prob_hi[dir]);

        if (beyond_lo || beyond_hi)
        {
            parity *= BoundaryConditions::get_state_var_parity(
                comp, dir); // Is there a boundary conditions function for
                            // derived vars? TODO!!!
            // invert parity for first derivatives
            if (deriv[dir] == 1)
                parity *= -1;
        }
    }
    return parity;
}

// a function to reflect a particle back into the valid domain, when symmetry
// BCs are used
amrex::Real ParticleInterpolators::reflect_particle(amrex::Real x,
                                                    amrex::Real lo,
                                                    amrex::Real hi,
                                                    bool lo_reflect,
                                                    bool hi_reflect) const
{
    // enforce a new particle position if needed
    amrex::Real xl = x;
    if (lo_reflect && xl < lo)
    {
        xl = lo + (lo - xl); // reflect across lo
        amrex::Print() << "Particle position " << x << " reflected to " << xl
                       << " across the low boundary.\n";
    }
    if (hi_reflect && xl > hi)
    {
        xl = hi - (xl - hi); // reflect across hi
        amrex::Print() << "Particle position " << x << " reflected to " << xl
                       << " across the high boundary.\n";
    }

    return xl;
}

// allocate particles at the query points
void ParticleInterpolators::populate_from_query(
    const InterpolationQueryParticle &query)
{
    AMREX_ASSERT(m_initialized);

    // we populate particles with rank 0
    const bool amroot = (amrex::ParallelDescriptor::MyProc() == 0);

    const int lev = 0;
    const int n =
        static_cast<int>(query.m_num_points); // number of query points

    if (amroot)
    {
        // it does not matter on which level the particles are initialised so
        // long as we Redistribute() them later
        auto &ptile = this->DefineAndReturnParticleTile(lev, 0, 0);
        ptile.resize(n);
        auto particle_data = ptile.getParticleTileData();

        // get coords from query
        const double *x = query.m_coords[0];
        const double *y = query.m_coords[1];
#if AMREX_SPACEDIM == 3
        const double *z = query.m_coords[2];
#endif

        // loop over particles and place them at the required points
        amrex::ParallelFor(
            n,
            [=] AMREX_GPU_DEVICE(int ip)
            {
                auto &p = particle_data[ip]; 
                p.id()  = ip + 1;            // particle id starts from 1
                p.cpu() = 0;                 // CPU id, not used here

                // reflect into valid region for seeding
                const amrex::Real xr = reflect_particle(
                    static_cast<amrex::Real>(x[ip]), m_prob_lo[0], m_prob_hi[0],
                    m_lo_boundary_reflective[0], m_hi_boundary_reflective[0]);
                const amrex::Real yr = reflect_particle(
                    static_cast<amrex::Real>(y[ip]), m_prob_lo[1], m_prob_hi[1],
                    m_lo_boundary_reflective[1], m_hi_boundary_reflective[1]);
#if AMREX_SPACEDIM == 3
                const amrex::Real zr = reflect_particle(
                    static_cast<amrex::Real>(z[ip]), m_prob_lo[2], m_prob_hi[2],
                    m_lo_boundary_reflective[2], m_hi_boundary_reflective[2]);
#endif

                // Run a check on coords you are interpolating on
                std::array<double, AMREX_SPACEDIM> coords{xr, yr
#if AMREX_SPACEDIM == 3
                                                          ,
                                                          zr
#endif
                };
                check_domain<AMREX_SPACEDIM>(coords, 0);

                // set position
                p.pos(0) = xr;
                p.pos(1) = yr;
#if AMREX_SPACEDIM == 3
                p.pos(2) = zr;
#endif
                p.idata(0) = ip;

                for (int s = 0; s < m_ncomp; ++s)
                    particle_data.rdata(s)[ip] = 0.0; // for now set all to zero (this is where we will store the interpolated values)
            });

        amrex::Gpu::streamSynchronize();
    }

    amrex::ParallelDescriptor::Barrier(); // TODO! test if this makes a
                                          // difference.
}

// interpolate variables into SOA slots
void ParticleInterpolators::interpolate_to_particle()
{
    Redistribute();

    AMREX_ASSERT(m_initialized);
    AMREX_ALWAYS_ASSERT(m_ncomp >= 1 && m_ncomp <= AMREX_SPACEDIM);
    AMREX_ASSERT(m_start_comp >= 0);

    for (int lev = 0; lev <= m_gr_amr->finestLevel(); ++lev)
    {
        // if particles are only on the coarsest level, skip the looping over
        // levels
        if (this->NumberOfParticlesAtLevel(lev) == 0)
            continue;

        amrex::AmrLevel &level      = m_gr_amr->getLevel(lev);
        const amrex::Geometry &geom = level.Geom();
        amrex::MultiFab &state      = level.get_new_data(State_Type);

        AMREX_ASSERT(m_start_comp + m_ncomp <= state.nComp());

        amrex::IntVect nghost(AMREX_D_DECL(2, 2, 2));
        state.FillBoundary(m_start_comp, m_ncomp, nghost, geom.periodicity());

        const auto plo = geom.ProbLoArray();
        const auto dxi = geom.InvCellSizeArray();

        // loop over tiles and interpolate now
        for (ParIterType it(*this, lev); it.isValid(); ++it)
        {
            auto &ptile    = this->ParticlesAt(lev, it);
            auto ptd       = ptile.getParticleTileData();
            const int np   = it.numParticles();
            auto fab_array = state[it].const_array();

            amrex::ParallelFor(
                np,
                [=] AMREX_GPU_DEVICE(int ip)
                {
                    auto &sp = ptd[ip];

                    amrex::IntVect is_nodal = amrex::IntVect::TheZeroVector();
                    LagrangeInterpolator<5> interp; // 4th order interpolation
                    interp.compute_weights(sp, plo, dxi, is_nodal);

                    amrex::ParticleReal vals[AMREX_SPACEDIM];
                    interp.interpolate(&fab_array, vals, m_start_comp, m_ncomp);

                    // write results to SOA
                    for (int k = 0; k < m_ncomp; ++k)
                    {
                        ptd.rdata(k)[ip] = vals[k];
                    }
                });
        }

        // synchronize GPU streams to ensure all particles are updated
        amrex::Gpu::streamSynchronize();
    }
}

// mirror of AMRInterpolator::interp(); assembles all particle data and writes
// parity * value into the query out arrays
void ParticleInterpolators::interp(InterpolationQueryParticle &query)
{
    AMREX_ASSERT(m_initialized);

    // get total query points here
    const int npts = static_cast<int>(query.numPoints());

    // value_at_point[k][ip], where k in [0..m_ncomp-1]
    std::vector<std::vector<double>> value_at_point(
        m_ncomp, std::vector<double>(npts, 0.0));
    // a vector to mark which points have values
    std::vector<int> have(npts, 0);

    int local_particle_counter = 0;
    // gather all particle values
    for (int lev = 0; lev <= this->finestLevel(); ++lev)
    {
        for (ParIterType it(*this, lev); it.isValid(); ++it)
        {
            local_particle_counter += it.numParticles();

            const auto &ptile = this->ParticlesAt(lev, it);
            const auto aos    = ptile.GetArrayOfStructs();
            const auto *P     = aos.data();
            const auto soa = ptile.getConstParticleTileData(); // read-only SoA
            const int np   = it.numParticles();

            for (int i = 0; i < np; ++i)
            {
                auto &p     = P[i];
                const int q = p.idata(0); // get particle index
                for (int k = 0; k < m_ncomp; ++k)
                {
                    value_at_point[k][q] = static_cast<double>(soa.rdata(k)[i]);
                }
                have[q] = 1; // mark that we have a value for this point
            }
        }
    }

    amrex::AllPrint() << "rank " << amrex::ParallelDescriptor::MyProc()
                      << " holds " << local_particle_counter << " particles\n";

    // reduce across ranks
    for (int k = 0; k < m_ncomp; ++k)
    {
        amrex::ParallelDescriptor::ReduceRealSum(value_at_point[k].data(),
                                                 npts);
    }
    amrex::ParallelDescriptor::ReduceIntSum(have.data(), npts);

    if (amrex::ParallelDescriptor::IOProcessor())
    {
        amrex::AllPrint() << "The IO rank is = "
                          << amrex::ParallelDescriptor::IOProcessorNumber()
                          << "\n";
        // loop for each component and apply parity assumptions now
        for (auto deriv_it = query.compsBegin(); deriv_it != query.compsEnd();
             ++deriv_it)
        {
            using comps_t =
                std::vector<typename InterpolationQueryParticle::out_t>;
            comps_t &comps = deriv_it->second;

            for (auto &entry : comps)
            {
                const int comp = std::get<0>(entry);
                double *out = std::get<1>(entry); // this is where interpolated
                                                  // values are written into
                const Derivative &dkey = deriv_it->first;

                const int k =
                    comp -
                    m_start_comp; // reindex the variable component from 0
                AMREX_ALWAYS_ASSERT(k >= 0 && k < m_ncomp);

                for (int ip = 0; ip < npts; ++ip)
                {
                    if (!have[ip])
                    {
                        amrex::Abort("interp(): no data for query point " +
                                     std::to_string(ip));
                    }

                    const int parity =
                        get_state_var_parity(comp, ip, query, dkey);
                    const double v = have[ip] ? value_at_point[k][ip] : 0.0;
                    out[ip]        = parity * v;
                }
            }
        }
    }
}

// A function to check whether the query point is inside the physical domain
template <int dim>
void ParticleInterpolators::check_domain(const std::array<double, dim> &x,
                                         int guard_cells)
{
    for (int d = 0; d < dim; ++d)
    {
        const double lo_g = m_prob_lo[d] - guard_cells * m_dx[d];
        const double hi_g = m_prob_hi[d] + guard_cells * m_dx[d];
        if (x[d] < lo_g || x[d] > hi_g)
        {
            std::string msg = "ParticleInterpolators::check_domain() Oi oi oi! "
                              "You are trying to access the point at x[" +
                              std::to_string(d) +
                              "] = " + std::to_string(x[d]) +
                              " and it lies outside of your domain.";

            amrex::Abort(msg);
        }
    }
}

// TODO: I have not tested the below yet!!

// Should I have a function to clear particles at a specific level?
void ParticleInterpolators::clear_level(int lev)
{
    for (ParIterType it(*this, lev); it.isValid(); ++it)
    {
        auto &ptile = this->ParticlesAt(lev, it);
        ptile.resize(0); // by resizing to 0 we clear the particles?
    }
    this->Redistribute(); // I am not sure about this :/
}

// Clear particles on all levels
void ParticleInterpolators::clear_all()
{
    for (int lev = 0; lev <= this->finestLevel(); ++lev)
    {
        clear_level(lev);
    }
}

#endif /* PARTICLEINTERPOLATORS_IMPL_HPP_ */
