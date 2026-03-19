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
    bool a_verbosity)
{
    m_tol = 1e-4;
    m_eta = 5;
    m_c   = 1.0;
    m_dt  = 0.01;

    amrex::Real r      = 1.15;
    amrex::Real min_dt = 1e-4;
    amrex::Real max_dt = 1e2;

    amrex::Real theta_old;
    amrex::Real theta_new;

    int n_iter = 0;

    // Create InterpolationQueryParticle
    // of n coordinates based on
    // centre and radius of guess
    this->generate_spherical_query();

    // Set up interpolator
    this->setup(gramr_ptr, a_bc_params, a_verbosity);

    // Populate so we can access the particle data
    this->populate_from_query();

    // Initialise h and v values for particles
    this->init_h_v();

    this->interp(query);

    theta_new = inf_norm(interp_vals);

    this->move_radial();

    while (theta_new > m_tol)
    {
        theta_old = theta_new;
        this->update_v();
        this->move_radial();
        this->interp(query);

        theta_new = inf_norm(interp_vals);

        // Update time step based on ratio of old to new theta
        // As we converge to theta = 0 we can increase the timestep
        m_dt = r * m_dt * theta_new / theta_old;

        // Ensure timestep doesn't grow too large or small
        m_dt = std::max(m_dt, min_dt);
        m_dt = std::min(m_dt, max_dt);

        n_iter++;
    }

    amrex::AllPrint() << "\n AHFinder converged with inf norm of theta = "
                      << theta_new << " in " << n_iter << " iterations\n";
}

template <int num_components>
void AHFinder<num_components>::generate_spherical_query()
{
    for (int j = 0; j < m_num_particles; ++j)
    {
        double phi   = j * 2. * M_PI / m_num_particles;
        double theta = j * M_PI / m_num_particles;

        interp_coords_x[j] = m_center[0] + cos(phi) * sin(theta);
        interp_coords_y[j] = m_center[1] + sin(phi) * sin(theta);
        interp_coords_z[j] = m_center[2] + cos(theta);
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
                                       sqrt(pow(p.pos(0) - m_center[0], 2) +
                                            pow(p.pos(1) - m_center[1], 2) +
                                            pow(p.pos(2) - m_center[2],
                                                2)); // Height from centre

                                   p.rdata(1) = 0.0; // Velocity
                               });
        }

        amrex::Gpu::streamSynchronize();
    }
}

template <int num_components> void AHFinder<num_components>::move_radial()
{
    for (int lev = 0; lev <= this->m_gramr_ptr->finestLevel(); ++lev)
    {

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
                    double center_direction;

                    auto &p = pstruct[ip];

                    amrex::GpuArray<double, AMREX_SPACEDIM> coords;

                    for (size_t i = 0; i < AMREX_SPACEDIM; i++)
                    {
                        // Get normalised direction to m_center
                        center_direction =
                            (p.pos(i) - m_center[i]) / p.rdata(0);

                        // Update position
                        p.pos(i) += m_dt * (center_direction *
                                            (p.rdata(1) - m_eta * p.rdata(0)));

                        coords[i] = p.pos(i);
                    }

                    this->check_domain(coords);

                    // Update h
                    p.rdata(0) = sqrt(pow(p.pos(0) - m_center[0], 2) +
                                      pow(p.pos(1) - m_center[1], 2) +
                                      pow(p.pos(2) - m_center[2], 2));

                    // Update query position
                    interp_coords_x[ip] = p.pos(0);
                    interp_coords_y[ip] = p.pos(1);
                    interp_coords_z[ip] = p.pos(2);
                });
        }

        amrex::Gpu::streamSynchronize();

        // Update query
        query.setCoords(0, interp_coords_x.data())
            .setCoords(1, interp_coords_y.data())
            .setCoords(2, interp_coords_z.data());
    }
}

template <int num_components> void AHFinder<num_components>::update_v()
{
    for (int lev = 0; lev <= this->m_gramr_ptr->finestLevel(); ++lev)
    {

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

                                   // Update v
                                   p.rdata(1) -=
                                       m_dt * pow(m_c, 2) * interp_vals[ip];
                               });
        }

        amrex::Gpu::streamSynchronize();
    }
}

template <int num_components>
double AHFinder<num_components>::inf_norm(std::vector<double> arr)
{
    double max_el = std::abs(arr[0]);
    for (auto &&i : arr)
    {
        if (std::abs(i) > max_el)
            max_el = std::abs(i);
    }

    return max_el;
}

#endif /* AHFINDER_IMPL_HPP_ */
