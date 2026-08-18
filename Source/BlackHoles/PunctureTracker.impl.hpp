/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef PUNCTURETRACKER_HPP_
#error "This file should only be included through PunctureTracker.hpp"
#endif

#ifndef PUNCTURETRACKER_IMPL_HPP_
#define PUNCTURETRACKER_IMPL_HPP_

// #include "AMReXParameters.hpp" // for writing data
#include "DimensionDefinitions.hpp"
#include "FilesystemTools.hpp"
#include "GRAMRLevel.hpp"
#include "SmallDataIO.hpp" // for writing data
#include "StateTypes.hpp"
#include "StateVariables.hpp"

// AMReX includes
#include <AMReX_AmrLevel.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_ParmParse.H>
#include <AMReX_TracerParticle_mod_K.H> // for linear_interpolation

void puncture_tracker_params_t::check_params()
{
    GRParmParse puncture_tracking_pp("puncture_tracking");
    GRParmParse pp;

    bool enabled = false; // default
    puncture_tracking_pp.queryAdd("enabled", enabled);
    if (!enabled)
    {
        return;
    }

    std::string output_path;
    pp.get("grteclyn.output_path", output_path);

    std::string pt_output_path = output_path + "/punctures_output";

    if (!FilesystemTools::directory_exists(pt_output_path))
    {
        FilesystemTools::mkdir_recursive(pt_output_path);
    }
    puncture_tracking_pp.add("output_path", pt_output_path);

    std::string filename = "punctures";
    puncture_tracking_pp.add("filename", filename);

    bool disable_writeout = false;
    puncture_tracking_pp.queryAdd("disable_writeout", disable_writeout);

    int max_level;
    pp.get("amr.max_level", max_level);
    int level = max_level;
    puncture_tracking_pp.queryAdd("level", level);
    if (level < 0 || level > max_level)
    {
        puncture_tracking_pp.warning(
            "level", "must be between 0 and max_level (inclusive)");
    }

    int writeout_level = 0;
    puncture_tracking_pp.queryAdd("writeout_level", writeout_level);

    std::array<double, AMREX_SPACEDIM> center{};
    pp.get("geometry.center", center);

    std::array<amrex::Real, AMREX_SPACEDIM * 2UL> initial_coords{
        center[0], center[1] - 1.0, center[2],
        center[0], center[1] + 1.0, center[2]};
    puncture_tracking_pp.queryAdd("initial_coords", initial_coords);
}

void puncture_tracker_params_t::fill_params()
{
    amrex::ParmParse puncture_tracking_pp("puncture_tracking");

    puncture_tracking_pp.get("output_path", output_path);
    puncture_tracking_pp.get("filename", filename);
    full_filename     = output_path + "/" + filename;
    checkpoint_subdir = filename;

    puncture_tracking_pp.get("disable_writeout", disable_writeout);

    puncture_tracking_pp.get("level", level);
    puncture_tracking_pp.get("writeout_level", writeout_level);

    puncture_tracking_pp.get("initial_coords", initial_coords);
}

//! Set up puncture tracker
template <unsigned int num_punctures>
void PunctureTracker<num_punctures>::initialize(GRAMR *a_gr_amr)
{
    m_params.fill_params();
    AMREX_ASSERT(a_gr_amr != nullptr);
    m_gr_amr = a_gr_amr;

    {
        // Disable particle tiling as we won't have many particles
        // TODO: Remove if we add more particles elsewhere
        amrex::ParmParse particles_pp("particles");
        particles_pp.add("do_tiling", 0);
    }

    m_initialized = true;
}

template <unsigned int num_punctures>
void PunctureTracker<num_punctures>::start_from_initial_punctures()
{
    AMREX_ASSERT(m_initialized);
    // must call set_puncture_coords for the initial punctures first
    AMREX_ASSERT(m_puncture_coords_set);

    // Define the particle container
    Define(dynamic_cast<amrex::ParGDBBase *>(m_gr_amr->GetParGDB()));

    // If it's first step, we use the initial puncture locations set above
    write_initial_punctures();

    // Add the initial puncture particles to the underlying
    // ParticleContainer
    set_initial_punctures_pc();

    m_started = true;
}

template <unsigned int num_punctures>
void PunctureTracker<num_punctures>::restart(
    const std::string &a_restart_chk_dir)
{
    AMREX_ASSERT(m_initialized);

    // Define the particle container
    Define(dynamic_cast<amrex::ParGDBBase *>(m_gr_amr->GetParGDB()));

    Restart(a_restart_chk_dir, m_params.checkpoint_subdir);

    m_started = true;

    m_restart_time = m_gr_amr->get_restart_time();

    // The above Restart function will only set the punctures in the underlying
    // ParticleContainer so let's update our own m_puncture_coords
    update_puncture_coords();
}

template <unsigned int num_punctures>
void PunctureTracker<num_punctures>::write_plotfile(const std::string &a_dir)
{
    AMREX_ASSERT(m_initialized);
    AMREX_ASSERT(m_started);

    std::string plotfile_subdir = "particles"; // this is what ParaView expects

    amrex::Vector<std::string> real_comp_names{AMREX_D_DECL(
        StateVariables::names[c_shift1], StateVariables::names[c_shift2],
        StateVariables::names[c_shift3])};

    amrex::Vector<std::string> int_comp_names({"puncture_index"});

    Redistribute();
    WritePlotFile(a_dir, plotfile_subdir, real_comp_names, int_comp_names);
}

template <unsigned int num_punctures>
void PunctureTracker<num_punctures>::checkpoint(const std::string &a_chk_dir)
{
    AMREX_ASSERT(m_initialized);
    AMREX_ASSERT(m_started);

    Redistribute();
    Checkpoint(a_chk_dir, m_params.checkpoint_subdir);
}

//! set and write initial puncture locations
template <unsigned int num_punctures>
void PunctureTracker<num_punctures>::set_initial_punctures_pc()
{
    AMREX_ASSERT(m_initialized);

    if (amrex::ParallelDescriptor::MyProc() != 0)
        return;

    // It doesn't matter where we put the puncture particles initially.
    // They will be redistributed later
    const int base_level = 0;
    {
        auto &particle_tile = DefineAndReturnParticleTile(base_level, 0, 0);
        particle_tile.resize(num_punctures);
        const auto &particle_tile_data = particle_tile.getParticleTileData();

        amrex::GpuArray<amrex::Real, num_puncture_coords> d_puncture_coords;
        std::copy(m_puncture_coords.begin(), m_puncture_coords.end(),
                  d_puncture_coords.begin());

        amrex::ParallelFor(
            num_punctures,
            [=] AMREX_GPU_DEVICE(int ipuncture)
            {
                FOR1 (idir)
                {
                    auto &puncture_particle = particle_tile_data[ipuncture];
                    puncture_particle.pos(idir) =
                        d_puncture_coords[linear_idx(ipuncture, idir)];
                    puncture_particle.id()     = ipuncture + 1;
                    puncture_particle.idata(0) = ipuncture + 1;
                    puncture_particle.cpu()    = 0;
                }
            });
        amrex::Gpu::streamSynchronize();
    }
}

template <unsigned int num_punctures>
void PunctureTracker<num_punctures>::set_puncture_coords(
    const amrex::Array<amrex::Real,
                       PunctureTracker<num_punctures>::num_puncture_coords>
        &a_puncture_coords)
{
    m_puncture_coords = a_puncture_coords;

    m_puncture_coords_set = true;
}

template <unsigned int num_punctures>
const amrex::Array<amrex::Real,
                   PunctureTracker<num_punctures>::num_puncture_coords> &
PunctureTracker<num_punctures>::get_puncture_coords() const
{
    AMREX_ASSERT(m_puncture_coords_set);
    return m_puncture_coords;
}

template <unsigned int num_punctures>
std::vector<amrex::Real>
PunctureTracker<num_punctures>::get_puncture_vector() const
{
    AMREX_ASSERT(m_initialized);
    AMREX_ASSERT(m_puncture_coords_set);

    std::vector<amrex::Real> puncture_coords_vector(num_puncture_coords);
    std::copy(m_puncture_coords.begin(), m_puncture_coords.end(),
              puncture_coords_vector.begin());

    return puncture_coords_vector;
}

template <unsigned int num_punctures>
void PunctureTracker<num_punctures>::write_initial_punctures() const
{
    AMREX_ASSERT(m_initialized);
    if (m_params.disable_writeout)
    {
        return;
    }
    // now the write out to a new file
    bool first_step = true;
    double dt       = 1.; // doesn't matter
    double time     = 0.;
    SmallDataIO punctures_file(m_params.full_filename, dt, time, m_restart_time,
                               SmallDataIO::APPEND, first_step);
    std::vector<std::string> header1_strings(
        static_cast<size_t>(num_puncture_coords));
    for (int ipuncture = 0; ipuncture < num_punctures; ipuncture++)
    {
        std::string idx = std::to_string(ipuncture + 1);
        header1_strings[AMREX_SPACEDIM * ipuncture + 0] = "x_" + idx;
        header1_strings[AMREX_SPACEDIM * ipuncture + 1] = "y_" + idx;
        header1_strings[AMREX_SPACEDIM * ipuncture + 2] = "z_" + idx;
    }
    punctures_file.write_header_line(header1_strings);

    // use a vector for the write out
    punctures_file.write_time_data_line(get_puncture_vector());
}

//! track the punctures and write out if requested
template <unsigned int num_punctures>
void PunctureTracker<num_punctures>::track(double a_time, double a_dt,
                                           const bool a_write_punctures)
{
    BL_PROFILE("PunctureTracker::track");
    AMREX_ASSERT(m_initialized);
    AMREX_ASSERT(m_started);

    // leave if this is called at t=0, we don't want to move the puncture yet
    {
        if (a_time == 0.)
            return;
    }

    // Redistribute punctures to the correct grid
    Redistribute();

    for (int ilevel = 0; ilevel <= m_gr_amr->finestLevel(); ilevel++)
    {
        if (this->NumberOfParticlesAtLevel(ilevel) == 0L)
        {
            continue;
        }
        amrex::AmrLevel &amr_level = m_gr_amr->getLevel(ilevel);

        const amrex::Geometry &geom  = amr_level.Geom();
        amrex::MultiFab &state_level = amr_level.get_new_data(state_index);

        // We should only need 1 ghost cell as we are doing linear interpolation
        amrex::IntVect ghosts_to_fill = amrex::IntVect::TheUnitVector();
        state_level.FillBoundary(c_shift1, GR_SPACEDIM, ghosts_to_fill,
                                 geom.periodicity());

        const auto problem_domain_lo = geom.ProbLoArray();
        const auto problem_domain_hi = geom.ProbHiArray();
        const auto dxi               = geom.InvCellSizeArray();

        // This code is almost identical to
        // TracerParticleContainer::AdvectWithUcc except we advect in the
        // opposite direction to the shift.
        for (int ipass = 0; ipass < 2; ipass++)
        {
            for (ParIterType punc_iter(*this, ilevel); punc_iter.isValid();
                 ++punc_iter)
            {
                ParticleTileType &punc_tile = ParticlesAt(ilevel, punc_iter);
                auto &punc_particles        = punc_tile.GetArrayOfStructs();
                auto *punc_particles_data   = punc_particles.data();
                int num_punc_tile           = punc_iter.numParticles();
                const auto &fab_array = state_level[punc_iter].const_array();

                amrex::ParallelFor(
                    num_punc_tile,
                    [=] AMREX_GPU_DEVICE(int ipunc)
                    {
                        auto &p = punc_particles_data[ipunc];
                        amrex::ParticleReal shift[AMREX_SPACEDIM];
                        amrex::IntVect is_nodal =
                            amrex::IntVect::TheZeroVector();
                        int num_arrays = 1;

                        amrex::linear_interpolate_to_particle(
                            p, problem_domain_lo, dxi, &fab_array, shift,
                            &is_nodal, c_shift1, GR_SPACEDIM, num_arrays);

                        if (ipass == 0)
                        {
                            FOR1 (idir)
                            {
                                p.rdata(idir) = p.pos(idir);
                                p.pos(idir) -= static_cast<amrex::ParticleReal>(
                                    0.5 * a_dt * shift[idir]);
                            }
                        }
                        else
                        {
                            FOR1 (idir)
                            {
                                p.pos(idir) = p.rdata(idir) -
                                              static_cast<amrex::ParticleReal>(
                                                  a_dt * shift[idir]);
                                p.rdata(idir) = shift[idir];
                            }
                        }

                        // make sure the particles don't leave the problem
                        // domain otherwise AMReX will mark them invalid
                        FOR1 (idir)
                        {
                            p.pos(idir) =
                                std::max(p.pos(idir), problem_domain_lo[idir]);
                            p.pos(idir) =
                                std::min(p.pos(idir), problem_domain_hi[idir]);
                        }
                    }); // amrex::ParallelFor
            } // punc_iter
        } // ipass
    } // ilevel

    // update m_puncture_coords with the updated locations of the puncture
    // particles
    update_puncture_coords();

    // write them out
    if (a_write_punctures && !m_params.disable_writeout)
    {
        bool first_step = false;
        SmallDataIO punctures_file(m_params.full_filename, a_dt, a_time,
                                   m_restart_time, SmallDataIO::APPEND,
                                   first_step);

        // use a vector for the write out
        punctures_file.remove_duplicate_time_data();
        punctures_file.write_time_data_line(get_puncture_vector());
    }
}

template <unsigned int num_punctures>
void PunctureTracker<num_punctures>::update_puncture_coords()
{
    BL_PROFILE("PunctureTracker::update_puncture_coords");
    AMREX_ASSERT(m_initialized);
    AMREX_ASSERT(m_started);

    // We will perform an MPI sum reduction to get the coords so set them to
    // zero by default
    m_puncture_coords.fill(0.0);

    amrex::Gpu::DeviceVector<amrex::ParticleReal> d_puncture_coords(
        num_puncture_coords, 0.0);

    auto d_puncture_coords_ptr = d_puncture_coords.data();

    for (int ilevel = 0; ilevel <= m_gr_amr->finestLevel(); ilevel++)
    {
        if (this->NumberOfParticlesAtLevel(ilevel) == 0L)
        {
            continue;
        }

        // Now if this proc has a puncture particle, we set its location
        for (ParIterType punc_iter(*this, ilevel); punc_iter.isValid();
             ++punc_iter)
        {
            auto &punc_particles      = punc_iter.GetArrayOfStructs();
            auto *punc_particles_data = punc_particles.data();
            int num_punc_tile         = punc_iter.numParticles();

            amrex::ParallelFor(
                num_punc_tile,
                [=] AMREX_GPU_DEVICE(int ipunc)
                {
                    auto &p      = punc_particles_data[ipunc];
                    int punc_idx = p.idata(0) - 1;
                    FOR1 (idir)
                    {
                        d_puncture_coords_ptr[linear_idx(punc_idx, idir)] +=
                            p.pos(idir);
                    }
                });
        }
    } // ilevel

    amrex::Gpu::copy(amrex::Gpu::deviceToHost, d_puncture_coords.begin(),
                     d_puncture_coords.end(), m_puncture_coords.data());

    // MPI sum over all ranks
    amrex::ParallelAllReduce::Sum(m_puncture_coords.data(), num_puncture_coords,
                                  amrex::ParallelContext::CommunicatorAll());
}

#endif
