#if !defined(AHFINDER_HPP_)
#error "This file should only be included through AHFinder.hpp"
#endif

#ifndef AHFINDER_IMPL_HPP_
#define AHFINDER_IMPL_HPP_

#include <AMReX_Array.H>
#include <AMReX_ParIter.H>
#include <AMReX_Particles.H>

#include "DefaultLevelFactory.hpp"
#include "Derivative.hpp"
#include "GRAMR.hpp"
#include "ParticleInterpolator.hpp"

template <int num_components>
void AHFinder<num_components>::init(
    GRAMR *gramr_ptr, const BoundaryConditions::params_t &a_bc_params,
    std::array<double, AMREX_SPACEDIM> center, bool a_verbosity)
{
    // Create InterpolationQueryParticle
    // of n coordinates based on
    // centre and radius of guess
    this->generate_spherical_query(center);

    // Set up interpolator
    this->setup(gramr_ptr, a_bc_params, a_verbosity);

    this->populate_from_query();

    // Initialise h and v values for particles
    this->init_h_v();

    this->move_radial();
}

template <int num_components>
void AHFinder<num_components>::generate_spherical_query(
    std::array<double, AMREX_SPACEDIM> center)
{

    for (int j = 0; j < m_num_particles; ++j)
    {
        double phi   = j * 2. * M_PI / m_num_particles;
        double theta = j * M_PI / m_num_particles;

        interp_coords_x[j] = center[0] + cos(phi) * sin(theta);
        interp_coords_y[j] = center[1] + sin(phi) * sin(theta);
        interp_coords_z[j] = center[2] + cos(theta);
    }

    query.setCoords(0, interp_coords_x.data())
        .setCoords(1, interp_coords_y.data())
        .setCoords(2, interp_coords_z.data())
        .addComp(0, interp_vals.data(), VariableType::state);

    this->m_query = &query;
}

template <int num_components> void AHFinder<num_components>::init_h_v()
{
    for (int lev = 0; lev <= this->m_gramr_ptr->finestLevel(); ++lev)
    {
        auto &level      = this->m_gramr_ptr->getLevel(lev);
        const auto &geom = level.Geom();

        if (this->NumberOfParticlesAtLevel(lev) == 0)
            continue;

        for (ParIterType par_iter(*this, lev); par_iter.isValid(); ++par_iter)
        {
            // Get AoS data for particles at this level
            auto &particle_tile   = this->ParticlesAt(lev, par_iter);
            auto &aos             = particle_tile.GetArrayOfStructs();
            ParticleType *pstruct = aos().dataPtr();

            amrex::ParallelFor(m_num_particles,
                               [=] AMREX_GPU_DEVICE(int ip)
                               {
                                   auto &p = pstruct[ip];

                                   p.rdata(0) =
                                       sqrt(pow(p.pos(0) - center[0], 2) +
                                            pow(p.pos(1) - center[1], 2) +
                                            pow(p.pos(2) - center[2],
                                                2)); // Height from centre
                                   p.rdata(1) = 1.0; // Velocity
                               });
        }

        amrex::Gpu::streamSynchronize();
    }
}

template <int num_components> void AHFinder<num_components>::move_radial()
{
    for (int lev = 0; lev <= this->m_gramr_ptr->finestLevel(); ++lev)
    {
        auto &level      = this->m_gramr_ptr->getLevel(lev);
        const auto &geom = level.Geom();

        if (this->NumberOfParticlesAtLevel(lev) == 0)
            continue;

        for (ParIterType par_iter(*this, lev); par_iter.isValid(); ++par_iter)
        {
            // Get AoS data for particles at this level
            auto &particle_tile   = this->ParticlesAt(lev, par_iter);
            auto &aos             = particle_tile.GetArrayOfStructs();
            ParticleType *pstruct = aos().dataPtr();

            amrex::ParallelFor(
                m_num_particles,
                [=] AMREX_GPU_DEVICE(int ip)
                {
                    amrex::Real center_direction;

                    auto &p = pstruct[ip];

                    for (size_t i = 0; i < AMREX_SPACEDIM; i++)
                    {
                        // Get direction to center
                        center_direction = p.pos(i) / p.rdata(0);

                        // Update position
                        p.pos(i) = p.pos(i) + center_direction * p.rdata(1);
                    }

                    // Update h
                    p.rdata(0) = sqrt(pow(p.pos(0) - center[0], 2) +
                                      pow(p.pos(1) - center[1], 2) +
                                      pow(p.pos(2) - center[2], 2));
                });
        }

        amrex::Gpu::streamSynchronize();
    }
}

#endif /* AHFINDER_IMPL_HPP_ */
