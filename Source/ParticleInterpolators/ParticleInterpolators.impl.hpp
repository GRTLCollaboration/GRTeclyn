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
#include "StateVariables.hpp"
#include "VariableType.hpp"

// amrex includes

#include <AMReX_AmrLevel.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>

template <int num_components>
ParticleInterpolators<num_components>::ParticleInterpolators(
    const BoundaryConditions::params_t &a_bc_params, int a_start_comp,
    bool a_verbosity)
    : m_gr_amr(nullptr), m_initialized(false), m_start_comp(a_start_comp),
      m_bc_params(a_bc_params), m_verbosity(a_verbosity)
{
}

// initialise everything and perform some sanity checks
template <int num_components>
void ParticleInterpolators<num_components>::set_gramr_ptr(GRAMR *gr_amr_ptr)
{
    // is GRAMR properly set?
    AMREX_ASSERT(gr_amr_ptr != nullptr);
    m_gr_amr = gr_amr_ptr;

    this->Define(dynamic_cast<amrex::ParGDBBase *>(m_gr_amr->GetParGDB()));
    m_initialized = true;

    AMREX_ALWAYS_ASSERT(num_components >= 1);
    AMREX_ALWAYS_ASSERT(m_start_comp >= 0);

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

template <int num_components>
void ParticleInterpolators<num_components>::set_derived_var_parity(int comp,
                                                                   BCParity p)
{
    AMREX_ALWAYS_ASSERT(comp >= 0 && comp < num_components);
    m_derived_bc_parity[comp] = p;
}

// a parity helper (the same way as it was defined in the AMRInterpolator)
template <int num_components>
int ParticleInterpolators<num_components>::get_state_var_parity(
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
ParticleInterpolators<num_components>::reflect_particle(amrex::Real x,
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
void ParticleInterpolators<num_components>::populate_from_query(
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

        // Run a check on coords you are interpolating on
        for (int i = 0; i < n; i++)
        {
            std::array<double, AMREX_SPACEDIM> coords{x[i], y[i]
#if AMREX_SPACEDIM == 3
                                                      ,
                                                      z[i]
#endif
            };
            check_domain(coords, 0);
        }

        // copy stuff as otherwise I get errors with lambdas when I use class
        // members inside the GPU loops; this is a bit annoying as I have to
        // copy a few things; any other alternative ways?
        const auto prob_lo    = m_prob_lo;
        const auto prob_hi    = m_prob_hi;
        const auto lo_reflect = m_lo_boundary_reflective;
        const auto hi_reflect = m_hi_boundary_reflective;
        const int ncomp       = num_components;

        // copy coords here
        const amrex::Real *x_p = nullptr, *y_p = nullptr, *z_p = nullptr;
        amrex::Gpu::DeviceVector<amrex::Real> Xd, Yd, Zd;

        // copy x,y,z to device vectors now
        amrex::Gpu::DeviceVector<double> X(n), Y(n);
        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, x, x + n, X.begin());
        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, y, y + n, Y.begin());

#if AMREX_SPACEDIM == 3
        amrex::Gpu::DeviceVector<double> Z(n);
        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, z, z + n, Z.begin());
#endif
        amrex::Gpu::streamSynchronize(); // ensure copies complete

        // Get raw device pointers to capture by value
        auto x_d = X.data();
        auto y_d = Y.data();
#if AMREX_SPACEDIM == 3
        auto z_d = Z.data();
#endif

        // loop over particles and place them at the required points
        amrex::ParallelFor(
            n,
            [=] AMREX_GPU_DEVICE(int ip)
            {
                auto &p = particle_data[ip];
                p.id()  = ip + 1; // particle id starts from 1
                p.cpu() = 0;      // CPU id, not used here

                // reflect into valid region for seeding
                const amrex::Real xr = reflect_particle(
                    static_cast<amrex::Real>(x_d[ip]), prob_lo[0], prob_hi[0],
                    lo_reflect[0], hi_reflect[0]);
                const amrex::Real yr = reflect_particle(
                    static_cast<amrex::Real>(y_d[ip]), prob_lo[1], prob_hi[1],
                    lo_reflect[1], hi_reflect[1]);
#if AMREX_SPACEDIM == 3
                const amrex::Real zr = reflect_particle(
                    static_cast<amrex::Real>(z_d[ip]), prob_lo[2], prob_hi[2],
                    lo_reflect[2], hi_reflect[2]);
#endif

                // set position
                p.pos(0) = xr;
                p.pos(1) = yr;
#if AMREX_SPACEDIM == 3
                p.pos(2) = zr;
#endif
                p.idata(0) = ip;

                for (int s = 0; s < ncomp; ++s)
                    particle_data.rdata(s)[ip] =
                        0.0; // for now set all to zero (this is where we will
                             // store the interpolated values)
            });

        amrex::Gpu::streamSynchronize();
    }

    // amrex::ParallelDescriptor::Barrier(); // TODO! test if this makes a
    // difference.
}

// interpolate variables into SOA slots
template <int num_components>
void ParticleInterpolators<num_components>::interpolate_to_particle()
{
    AMREX_ASSERT(m_initialized);

    ensure_redistributed();

    const int start_comp = m_start_comp;
    const int ncomp      = num_components;

    for (int lev = 0; lev <= m_gr_amr->finestLevel(); ++lev)
    {
        // if particles are only on the coarsest level, skip the looping over
        // levels
        if (this->NumberOfParticlesAtLevel(lev) == 0)
            continue;

        amrex::AmrLevel &level      = m_gr_amr->getLevel(lev);
        const amrex::Geometry &geom = level.Geom();
        amrex::MultiFab &state      = level.get_new_data(static_cast<int>(0));

        AMREX_ASSERT(start_comp + ncomp <= state.nComp());

        amrex::IntVect nghost(AMREX_D_DECL(2, 2, 2));
        state.FillBoundary(start_comp, ncomp, nghost, geom.periodicity());

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
                    interp.interpolate(&fab_array, vals, start_comp, ncomp);

                    // write results to SOA
                    for (int k = 0; k < ncomp; ++k)
                    {
                        ptd.rdata(k)[ip] = vals[k];
                    }
                });
        }

        // synchronize GPU streams to ensure all particles are updated
        amrex::Gpu::streamSynchronize();
    }

    m_particles_seeded  = true;
    m_need_redistribute = false;
}

// Interpolation for derived vars, takes in MultiFab and comps (unique
// components numbers), as e.g. we may want to interpolate several fields at
// once. However, as per current implementation, we can have only contiguous
// comps
template <int num_components>
void ParticleInterpolators<num_components>::
    interpolate_to_particle_from_derived_fields(
        const std::vector<const amrex::MultiFab *> &fields_by_lev)
{
    AMREX_ASSERT(m_initialized);

    const int nlevs = m_gr_amr->finestLevel() + 1;
    AMREX_ASSERT((int)fields_by_lev.size() == nlevs);

    ensure_redistributed();

    const int start_comp = m_start_comp;
    const int ncomp      = num_components;

    for (int lev = 0; lev <= m_gr_amr->finestLevel(); ++lev)
    {
        if (this->NumberOfParticlesAtLevel(lev) == 0)
            continue;

        auto &level               = m_gr_amr->getLevel(lev);
        const auto &geom          = level.Geom();
        const amrex::MultiFab &mf = *fields_by_lev[lev];

        AMREX_ASSERT(mf.nComp() >= start_comp + ncomp);

        // Fill boundaries
        amrex::IntVect nghost(AMREX_D_DECL(2, 2, 2));
        const_cast<amrex::MultiFab &>(mf).FillBoundary(
            start_comp, ncomp, nghost, geom.periodicity());

        const auto plo = geom.ProbLoArray();
        const auto dxi = geom.InvCellSizeArray();

        for (ParIterType it(*this, lev); it.isValid(); ++it)
        {
            auto arrs    = mf[it].const_array();
            auto &ptile  = this->ParticlesAt(lev, it);
            auto ptd     = ptile.getParticleTileData();
            const int np = it.numRealParticles();

            amrex::ParallelFor(
                np,
                [=] AMREX_GPU_DEVICE(int ip)
                {
                    auto &sp = ptd[ip];

                    amrex::IntVect is_nodal = amrex::IntVect::TheZeroVector();
                    LagrangeInterpolator<5> interp;
                    interp.compute_weights(sp, plo, dxi, is_nodal);

                    amrex::ParticleReal vals[ncomp];
                    interp.interpolate(&arrs, vals, start_comp, ncomp);
                    for (int k = 0; k < ncomp; ++k)
                    {
                        ptd.rdata(k)[ip] = vals[k];
                    }
                });

            amrex::Gpu::streamSynchronize();
        }
    }

    m_particles_seeded  = true;
    m_need_redistribute = false;
}

// mirror of AMRInterpolator::interp(); assembles all particle data and writes
// parity * value into the query out arrays
// Throws an error now if use requests interpolation of a derived variable and
// applies reflective BCs (needs fixing later!!) This is currently done on CPUs;
// I am not sure if there is an advantage of having any of this on GPUs
template <int num_components>
void ParticleInterpolators<num_components>::interp(
    InterpolationQueryParticle &query, VariableType variable_type)
{
    AMREX_ASSERT(m_initialized);

    // if (variable_type == VariableType::derived ||
    //     variable_type == VariableType::state)
    // {
    //     for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
    //     {
    //         if (m_lo_boundary_reflective[dir] ||
    //         m_hi_boundary_reflective[dir])
    //         {
    //             amrex::Abort("Help!!! reflective BCs are implemented, but not
    //             "
    //                          "yet tested. Likely to give absolute
    //                          gibberish!");
    //         }
    //     }
    // }

    // get total query points here
    const int npts = static_cast<int>(query.numPoints());

    // value_at_point[k][ip], where k in [0..num_components-1]
    std::vector<std::vector<amrex::Real>> value_at_point(
        num_components, std::vector<amrex::Real>(npts, 0.0));
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
            const auto soa    = ptile.getConstParticleTileData();
            const int np      = it.numParticles();

            amrex::Gpu::HostVector<ParticleType> hp(np);
            amrex::Gpu::copyAsync(amrex::Gpu::deviceToHost, aos.data(),
                                  aos.data() + np, hp.begin());

            std::vector<amrex::Gpu::HostVector<amrex::ParticleReal>> hr(
                num_components);

            for (int k = 0; k < num_components; ++k)
            {
                hr[k].resize(np);
                auto dptr = ptile.GetStructOfArrays().GetRealData(k).data();
                amrex::Gpu::copyAsync(amrex::Gpu::deviceToHost, dptr, dptr + np,
                                      hr[k].begin());
            }
            amrex::Gpu::streamSynchronize();

            for (int i = 0; i < np; ++i)
            {
                const int pid = hp[i].id();
                const int q   = static_cast<int>(pid) - 1; // get particle index
                if (q < 0 || q >= npts)
                {
                    amrex::Abort("interp(): particle id out of range");
                }
                for (int k = 0; k < num_components; ++k)
                {
                    value_at_point[k][q] = static_cast<double>(hr[k][i]);
                }
                have[q] = 1; // mark that we have a value for this point
            }
        }
    }
    if (m_verbosity)
    {
        amrex::AllPrint() << "rank " << amrex::ParallelDescriptor::MyProc()
                          << " holds " << local_particle_counter
                          << " particles\n";
    }

    // reduce across ranks
    for (int k = 0; k < num_components; ++k)
    {
        amrex::ParallelDescriptor::ReduceRealSum(value_at_point[k].data(),
                                                 npts);
    }
    amrex::ParallelDescriptor::ReduceIntSum(have.data(), npts);

    if (amrex::ParallelDescriptor::IOProcessor())
    {
        if (m_verbosity)
        {
            amrex::Print() << "m_comps size of the query is = "
                           << query.numComps() << "\n";
        }

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
                    m_start_comp; // reindex the variable component from 0;
                                  // works only for contiguous components
                AMREX_ALWAYS_ASSERT(k >= 0 && k < num_components);

                for (int ip = 0; ip < npts; ++ip)
                {
                    if (!have[ip])
                    {
                        amrex::Abort("interp(): no data for query point " +
                                     std::to_string(ip));
                    }

                    int parity = get_state_var_parity(comp, ip, query, dkey,
                                                      variable_type);
                    // amrex::Print() << "Point " << ip << " comp " << comp
                    //                << " value " << value_at_point[k][ip]
                    //                << " has parity " << parity << "\n";

                    const double v = have[ip] ? value_at_point[k][ip] : 0.0;
                    out[ip]        = parity * v;
                }
            }
        }
    }
}

// A function to check whether the query point is inside the physical domain
template <int num_components>
void ParticleInterpolators<num_components>::check_domain(
    std::array<double, AMREX_SPACEDIM> &x, int guard_cells) const
{
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
    {

        const double lo_g = m_prob_lo[d] - guard_cells * m_dx[d];
        const double hi_g = m_prob_hi[d] + guard_cells * m_dx[d];

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
        if (x[d] < m_prob_lo[d] && !m_lo_boundary_reflective[d] && x[d] < lo_g)
            error = true;
        if (x[d] > m_prob_hi[d] && !m_hi_boundary_reflective[d] && x[d] > hi_g)
            error = true;

        // still outside physical domain, even after applying reflective
        // conditions
        if (xr < lo_g || xr > hi_g)
            error = true;

        if (error)
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

// Ensure that particles are redistributed if needed
template <int num_components>
void ParticleInterpolators<num_components>::ensure_redistributed()
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
void ParticleInterpolators<num_components>::force_redistribute(bool flag)
{
    m_need_redistribute = flag;
}

// TODO: I have not tested the below yet!!

// Should I have a function to clear particles at a specific level?
template <int num_components>
void ParticleInterpolators<num_components>::clear_level(int lev)
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
void ParticleInterpolators<num_components>::clear_all()
{
    for (int lev = 0; lev <= this->finestLevel(); ++lev)
    {
        clear_level(lev);
    }
}

#endif /* PARTICLEINTERPOLATORS_IMPL_HPP_ */
