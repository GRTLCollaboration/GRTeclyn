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

    // Add comps for h and v
    this->AddRealComp(true);
    this->AddRealComp(true);

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
    // Generate particle coordinates such that they are laid out
    // in rings on a sphere

    // Aim for ~twice as many longitudes as latitudes
    m_n_rings  = 1;
    int target = std::max(1, (int)std::round(std::sqrt(m_num_particles / 2.0)));
    for (int n = 1; n <= m_num_particles; ++n)
    {
        if (m_num_particles % n == 0 &&
            std::abs(n - target) < std::abs(m_n_rings - target))
        {
            m_n_rings = n;
        }
    }
    m_ring_size = m_num_particles / m_n_rings;

    for (int i = 0; i < m_n_rings; ++i)
    {
        // Offset theta so we don't get points on the poles
        double theta = (i + 0.5) * M_PI / m_n_rings;
        for (int j = 0; j < m_ring_size; ++j)
        {
            double phi = j * 2. * M_PI / m_ring_size;
            int idx    = i * m_ring_size + j;

            interp_coords_x[idx] = m_center[0] + cos(phi) * sin(theta);
            interp_coords_y[idx] = m_center[1] + sin(phi) * sin(theta);
            interp_coords_z[idx] = m_center[2] + cos(theta);
        }
    }

    query.setCoords(0, interp_coords_x.data())
        .setCoords(1, interp_coords_y.data())
        .setCoords(2, interp_coords_z.data())
        .addComp(0, interp_vals.data(), VariableType::state);

    this->m_query = &query;
}

template <int num_components>
std::array<int, 4> AHFinder<num_components>::neighbours(int j) const
{
    // Get the 4 neighbours of a particle (north, south, east and west)

    int ring_num = j / m_ring_size;
    int ring_pos = j % m_ring_size;

    int north, south, east, west;

    // North/south neighbours are next ring above/below
    // In the case we are at a pole, take the particle opposite
    // on the same ring
    if (ring_num > 0)
    {
        north = (ring_num - 1) * m_ring_size + ring_pos;
    }
    else
    {
        north = (ring_num + (m_ring_size / 2)) % m_ring_size;
    }

    if (ring_num < m_n_rings - 1)
    {
        south = (ring_num + 1) * m_ring_size + ring_pos;
    }
    else
    {
        south = (ring_num + (m_ring_size / 2)) % m_ring_size;
    }

    // East/west neighbours are adjacent in the same ring
    // Modulo to wrap around
    east = (ring_num * m_ring_size + (ring_pos + 1)) % m_ring_size;
    west = ring_num * m_ring_size + (ring_pos - 1 + m_ring_size) % m_ring_size;

    return {north, south, east, west};
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
            auto &soa             = particle_tile.GetStructOfArrays();
            auto &aos             = particle_tile.GetArrayOfStructs();
            ParticleType *pstruct = aos().dataPtr();

            double *h_ptr = soa.GetRealData(m_h_idx).dataPtr();
            double *v_ptr = soa.GetRealData(m_v_idx).dataPtr();

            amrex::ParallelFor(m_num_particles,
                               [=] AMREX_GPU_DEVICE(int ip)
                               {
                                   auto &p = pstruct[ip];

                                   h_ptr[ip] =
                                       sqrt(pow(p.pos(0) - m_center[0], 2) +
                                            pow(p.pos(1) - m_center[1], 2) +
                                            pow(p.pos(2) - m_center[2],
                                                2)); // Height from centre

                                   v_ptr[ip] = 0.0; // Velocity
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
            auto &soa             = particle_tile.GetStructOfArrays();
            auto &aos             = particle_tile.GetArrayOfStructs();
            ParticleType *pstruct = aos().dataPtr();

            double *h_ptr = soa.GetRealData(m_h_idx).dataPtr();
            double *v_ptr = soa.GetRealData(m_v_idx).dataPtr();

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
                        center_direction = (p.pos(i) - m_center[i]) / h_ptr[ip];

                        // Update position
                        p.pos(i) += m_dt * (center_direction *
                                            (v_ptr[ip] - m_eta * h_ptr[ip]));

                        coords[i] = p.pos(i);
                    }

                    this->check_domain(coords);

                    // Update h
                    h_ptr[ip] = sqrt(pow(p.pos(0) - m_center[0], 2) +
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
            auto &particle_tile = this->ParticlesAt(lev, par_iter);
            auto &soa           = particle_tile.GetStructOfArrays();
            double *v_ptr       = soa.GetRealData(m_v_idx).dataPtr();

            amrex::ParallelFor(m_num_particles,
                               [=] AMREX_GPU_DEVICE(int ip)
                               {
                                   // Update v
                                   v_ptr[ip] -=
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
