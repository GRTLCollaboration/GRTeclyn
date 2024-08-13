/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#include "PunctureTracker.hpp"
// #include "AMReXParameters.hpp" // for writing data
#include "DimensionDefinitions.hpp"
#include "FilesystemTools.hpp"
#include "GRAMRLevel.hpp"
#include "GRParmParse.hpp"
#include "SmallDataIO.hpp"   // for writing data
#include "UserVariables.hpp" // for writing data

// AMReX includes
#include <AMReX_AmrParGDB.H>
#include <AMReX_TracerParticle_mod_K.H> // for linear_interpolation

//! Set up puncture tracker
void PunctureTracker::initial_setup(
    const amrex::Vector<amrex::RealVect> &initial_puncture_coords,
    GRAMR *a_gr_amr, const std::string &a_filename,
    const std::string &a_output_path, const int a_update_level)
{
    if (!FilesystemTools::directory_exists(a_output_path))
    {
        FilesystemTools::mkdir_recursive(a_output_path);
    }

    m_punctures_filename = a_output_path + a_filename;
    m_checkpoint_subdir  = a_filename;

    AMREX_ASSERT(a_gr_amr != nullptr);
    m_gr_amr = a_gr_amr;

    m_num_punctures     = static_cast<int>(initial_puncture_coords.size());
    m_puncture_coords   = initial_puncture_coords;
    m_puncture_proc_ids = amrex::Vector<int>(m_num_punctures, 0);

    m_update_level = a_update_level;

    {
        // Disable particle tiling as we won't have many particles
        // TODO: Remove if we add more particles elsewhere
        GRParmParse pp("particles");
        pp.add("do_tiling", 0);
    }
}

void PunctureTracker::restart(int a_coarse_step)
{
    Define(dynamic_cast<amrex::ParGDBBase *>(m_gr_amr->GetParGDB()));
    if (a_coarse_step == 0)
    {
        // if it is the first timestep, use the param values
        set_initial_punctures();
    }
    else
    {
        std::string restart_checkpoint{};
        GRParmParse pp("amr");
        pp.get("restart", restart_checkpoint);

        Restart(restart_checkpoint, m_checkpoint_subdir);
    }
}

void PunctureTracker::checkpoint(const std::string &a_chk_dir)
{
    redistribute();
    Checkpoint(a_chk_dir, m_checkpoint_subdir);
}

//! set and write initial puncture locations
void PunctureTracker::set_initial_punctures()
{
    AMREX_ASSERT(m_puncture_coords.size() > 0); // sanity check

    // now the write out to a new file
    bool first_step     = true;
    double dt           = 1.; // doesn't matter
    double time         = 0.;
    double restart_time = 0.;
    SmallDataIO punctures_file(m_punctures_filename, dt, time, restart_time,
                               SmallDataIO::APPEND, first_step);
    std::vector<std::string> header1_strings(
        static_cast<size_t>(AMREX_SPACEDIM * m_num_punctures));
    for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
    {
        std::string idx = std::to_string(ipuncture + 1);
        header1_strings[AMREX_SPACEDIM * ipuncture + 0] = "x_" + idx;
        header1_strings[AMREX_SPACEDIM * ipuncture + 1] = "y_" + idx;
        header1_strings[AMREX_SPACEDIM * ipuncture + 2] = "z_" + idx;
    }
    punctures_file.write_header_line(header1_strings);

    // use a vector for the write out
    punctures_file.write_time_data_line(get_puncture_vector());

    if (amrex::ParallelDescriptor::MyProc() != 0)
        return;

    // It doesn't matter where we put the puncture particles initially.
    // They will be redistributed later
    const int base_level = 0;
    {
        auto &particle_tile = DefineAndReturnParticleTile(base_level, 0, 0);

        particle_tile.resize(m_num_punctures);

        const auto &particle_tile_data = particle_tile.getParticleTileData();

        for (int ipuncture = 0; ipuncture < m_num_punctures; ++ipuncture)
        {
            // Maybe we can get away with doing this on the host
            FOR1 (idir)
            {
                auto &puncture_particle = particle_tile_data[ipuncture];
                puncture_particle.pos(idir) =
                    m_puncture_coords[ipuncture][idir];
                puncture_particle.id()  = ipuncture + 1;
                puncture_particle.cpu() = 0;
            }
        }
    }
}

void PunctureTracker::redistribute()
{
    // First do AMReX's particle redistribute
    Redistribute();

    // Now figure out which process has each puncture
    amrex::Vector<int> local_proc_has_punctures(m_num_punctures, 0);

    for (int ilevel = 0; ilevel <= m_gr_amr->finestLevel(); ilevel++)
    {
        if (this->NumberOfParticlesAtLevel(ilevel) == 0L)
        {
            continue;
        }
        for (ParIterType punc_iter(*this, ilevel); punc_iter.isValid();
             ++punc_iter)
        {
            auto &punc_particles = punc_iter.GetArrayOfStructs();
            int num_punc_tile    = punc_iter.numParticles();

            for (int ipunc = 0; ipunc < num_punc_tile; ipunc++)
            {
                ParticleType &p                     = punc_particles[ipunc];
                int ipuncture                       = p.id() - 1;
                local_proc_has_punctures[ipuncture] = 1;
            }
        }
    }

    amrex::Vector<int> global_proc_has_punctures(
        m_num_punctures * amrex::ParallelDescriptor::NProcs(), 0);

    // Communicate whether I have a puncture to all processes
    amrex::ParallelAllGather::AllGather(
        local_proc_has_punctures.dataPtr(), m_num_punctures,
        global_proc_has_punctures.dataPtr(),
        amrex::ParallelDescriptor::Communicator());

    // Keep track of the total number of punctures across all processes
    int puncture_count = 0;

    for (int iproc = 0; iproc < amrex::ParallelDescriptor::NProcs(); iproc++)
    {
        for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
        {
            if (global_proc_has_punctures[m_num_punctures * iproc +
                                          ipuncture] == 1)
            {
                m_puncture_proc_ids[ipuncture] = iproc;
                puncture_count++;
            }
        }
    }
    AMREX_ALWAYS_ASSERT(puncture_count == m_num_punctures);
}

//! Execute the tracking and write out
void PunctureTracker::execute_tracking(double a_time, double a_restart_time,
                                       double a_dt, const bool write_punctures)
{
    BL_PROFILE("PunctureTracker::execute_tracking");
    // leave if this is called at t=0, we don't want to move the puncture yet
    {
        if (a_time == 0.)
            return;
    }

    AMREX_ASSERT(static_cast<int>(m_puncture_coords.size()) ==
                 m_num_punctures); // sanity check

    // Redistribute punctures to the correct grid
    redistribute();

    for (int ilevel = 0; ilevel <= m_gr_amr->finestLevel(); ilevel++)
    {
        if (this->NumberOfParticlesAtLevel(ilevel) == 0L)
        {
            continue;
        }
        amrex::AmrLevel &amr_level = m_gr_amr->getLevel(ilevel);

        const amrex::Geometry &geom  = amr_level.Geom();
        amrex::MultiFab &state_level = amr_level.get_new_data(State_Type);

        amrex::MultiFab shift_mf(state_level, amrex::make_alias, c_shift1,
                                 AMREX_SPACEDIM);

        // We should only need 1 ghost cell as we are doing linear interpolation
        amrex::IntVect ghosts_to_fill = amrex::IntVect::TheUnitVector();
        shift_mf.FillBoundary(ghosts_to_fill, geom.periodicity());

        const auto problem_domain_lo = geom.ProbLoArray();
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
                // const auto &fab_array = state_level[punc_iter].const_array();
                const auto &shift_array = shift_mf[punc_iter].const_array();

                amrex::ParallelFor(
                    num_punc_tile,
                    [=] AMREX_GPU_DEVICE(int ipunc)
                    {
                        auto &p = punc_particles_data[ipunc];
                        amrex::ParticleReal shift[AMREX_SPACEDIM];
                        amrex::IntVect is_nodal =
                            amrex::IntVect::TheZeroVector();
                        int num_arrays = 1;

                        // amrex::linear_interpolate_to_particle(
                        //     p, problem_domain_lo, dxi, &fab_array, shift,
                        //     &is_nodal, c_shift1, AMREX_SPACEDIM, num_arrays);
                        cic_interpolate(p, problem_domain_lo, dxi, shift_array,
                                        shift);

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
                    }); // amrex::ParallelFor
            }           // punc_iter
        }               // ipass

        for (ParIterType punc_iter(*this, ilevel); punc_iter.isValid();
             ++punc_iter)
        {
            auto &punc_particles = punc_iter.GetArrayOfStructs();
            int num_punc_tile    = punc_iter.numParticles();

            for (int ipunc = 0; ipunc < num_punc_tile; ipunc++)
            {
                ParticleType &p = punc_particles[ipunc];

                m_puncture_coords[p.id() - 1] = p.pos();
            }
        }
    } // ilevel

    // broadcast the locations to all ranks
    for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
    {
        // If there are lots of punctures, we should probably use non-blocking
        // MPI calls which are currently not wrapped by AMReX
        amrex::ParallelDescriptor::Bcast(m_puncture_coords[ipuncture].dataPtr(),
                                         AMREX_SPACEDIM,
                                         m_puncture_proc_ids[ipuncture]);
    }

    // print them out
    if (write_punctures)
    {
        bool first_step = false;
        SmallDataIO punctures_file(m_punctures_filename, a_dt, a_time,
                                   a_restart_time, SmallDataIO::APPEND,
                                   first_step);

        // use a vector for the write out
        punctures_file.write_time_data_line(get_puncture_vector());
    }
}

//! get a vector of the puncture coords - used for write out
std::vector<double> PunctureTracker::get_puncture_vector() const
{
    std::vector<double> puncture_vector;
    puncture_vector.resize(m_num_punctures * AMREX_SPACEDIM); // NOLINT
    for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
    {
        puncture_vector[ipuncture * AMREX_SPACEDIM + 0] =
            m_puncture_coords[ipuncture][0];
        puncture_vector[ipuncture * AMREX_SPACEDIM + 1] =
            m_puncture_coords[ipuncture][1];
        puncture_vector[ipuncture * AMREX_SPACEDIM + 2] =
            m_puncture_coords[ipuncture][2];
    }
    return puncture_vector;
}
