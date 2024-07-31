/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#include "PunctureTracker.hpp"
// #include "AMReXParameters.hpp" // for writing data
#include "DimensionDefinitions.hpp"
#include "GRParmParse.hpp"
#include "SmallDataIO.hpp"   // for writing data
#include "UserVariables.hpp" // for writing data

// AMReX includes
#include <AMReX_TracerParticle_mod_K.H>

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

    AMREX_ASSERT(a_amr != nullptr);
    m_amr = a_amr;

    // first set the puncture data
    m_puncture_coords = initial_puncture_coords;
    m_num_punctures   = static_cast<int>(m_puncture_coords.size());

    m_update_level = a_update_level;

    {
        // Disable particle tiling as we won't have many particles
        // TODO: Remove if we add interpolation particles elsewhere
        GRParmParse pp("particles");
        pp.add("do_tiling", 0);
    }

    s_particle_container =
        std::make_unique<amrex::ParticleContainerPureSoA<0, 0>>(
            dynamic_cast<amrex::ParGDBBase *>(m_gr_amr->GetParGDB()));

    // It doesn't matter where we put the puncture particles initially.
    // They will be redistributed later
    const int base_level = 0;
    for (amrex::MFIter mfi = s_particle_container->MakeMFIter(base_level);
         mfi.isValid(); ++mfi)
    {
        auto &particle_tile = s_particle_container->DefineAndReturnParticleTile(
            base_level, mfi.index(), mfi.LocalTileIndex());
        if (mfi.index() != 0 || mfi.LocalTileIndex() != 0)
            continue;

        particle_tile.resize(m_num_punctures);

        const auto &particle_tile_data = particle_tile.getParticleTileData();

        amrex::ParallelFor(m_num_punctures,
                           [=] AMREX_GPU_DEVICE(int ipuncture)
                           {
                               FOR1 (idim)
                               {

                                   particle_tile_data.pos(idim, ipuncture) =
                                   // TODO: get the location on the GPUs
                                   // m_puncture_coords[ipuncture][idim];
                                   // TODO: initialize the real data (shift) to
                                   // 0.0
                               }
                           });
    }
}

void PunctureTracker::restart(int a_coarse_step)
{

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

        s_particle_container->Restart(restart_checkpoint, m_checkpoint_subdir);
    }
}

void PunctureTracker::checkpoint(const std::string &a_chk_dir)
{
    s_particle_container->Checkpoint(a_chk_dir, m_checkpoint_subdir);
}

//! set and write initial puncture locations
void PunctureTracker::set_initial_punctures()
{
    AMREX_ASSERT(m_puncture_coords.size() > 0); // sanity check

    m_puncture_shift.resize(m_num_punctures);
    for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
    {
        // assume initial shift is always zero
        FOR (i)
        {
            m_puncture_shift[ipuncture][i] = 0.0;
        }
    }

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
    s_particle_container->Redistribute();

    for (int ilevel = 0; ilevel < m_gr_amr->finestLevel(); ilevel++)
    {
        if (s_particle_container->NumberOfParticlesAtLevel(ilevel) == 0L)
        {
            continue;
        }
        const AmrLevel &amr_level = m_gr_amr->getLevel(ilevel);

        amrex::MultiFab &state_level = amr_level.get_new_data(State_Type);
        state_level.FillBoundary(c_shift1, GR_SPACEDIM);

        const Geometry &geom         = amr_level.Geom();
        const auto problem_domain_lo = geom.ProbLo();
        const auto dxi               = geom.InvCellSizeArray();

        // This code is almost identical to
        // TracerParticleContainer::AdvectWithUcc except we advect in the
        // opposite direction to the shift.
        for (int ipass = 0; ipass < 2; ipass++)
        {
            for (PunctureIter punc_iter(s_particle_container, ilevel);
                 punc_iter.isValid(); ++punc_iter)
            {
                auto &punc_tile =
                    s_particle_container->ParticlesAt(ilevel, punc_iter);
                int num_punc_tile     = punc_tile.numParticles();
                const auto &fab_array = state_level.const_arrays()[punc_iter];

                amrex::ParallelFor(
                    num_punc_tile,
                    [=] AMREX_GPU_DEVICE(int ipunc)
                    {
                        amrex::PunctureParticleType &p = punc_tile[ipunc];
                        amrex::ParticleReal shift[AMREX_SPACEDIM];
                        amrex::IntVect is_nodal =
                            amrex::IntVect::TheZeroVector();
                        int num_arrays = 1;

                        amrex::linear_interpolate_to_particle(
                            p, problem_domain_lo, dxi, &fab_array, shift,
                            is_nodal, c_shift1, AMREX_SPACEIM, num_arrays);

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
                                p.pos(idir) =
                                    p.rdata(idir) - static_cast<ParticleReal>(
                                                        a_dt * shift[idir]);
                                p.rdata(idir) = shift[idir];
                            }
                        }
                    }); // amrex::ParallelFor

                // xxxxx copy the puncture location back to host
            } // punc_iter
        }     // ipass
        // TODO: copy the puncture location back to host
        // TODO: broadcast the location to all ranks
    }

    // for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
    // {
    //     FOR (i)
    //     {
    //         m_puncture_coords[ipuncture][i] +=
    //             -0.5 * a_dt *
    //             (m_puncture_shift[ipuncture][i] + old_shift[ipuncture][i]);
    //     }
    // }

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

//! Use the interpolator to get the value of the shift at
//! given coords
void PunctureTracker::interp_shift()
{
#if 0
//xxxxx
    BL_PROFILE("PunctureTracker::interp_shift");
    // resize the vector to the number of punctures
    m_puncture_shift.resize(m_num_punctures);

    // refresh interpolator
    bool fill_ghosts = false;
    m_interpolator->refresh(fill_ghosts);
    // only fill the ghosts we need
    m_interpolator->fill_multilevel_ghosts(
        VariableType::evolution, Interval(c_shift1, c_shift3), m_min_level);

    // set up shift and coordinate holders
    std::vector<double> interp_shift1(m_num_punctures);
    std::vector<double> interp_shift2(m_num_punctures);
    std::vector<double> interp_shift3(m_num_punctures);
    std::vector<double> interp_x(m_num_punctures);
    std::vector<double> interp_y(m_num_punctures);
    std::vector<double> interp_z(m_num_punctures);

    for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
    {
        interp_x[ipuncture] = m_puncture_coords[ipuncture][0];
        interp_y[ipuncture] = m_puncture_coords[ipuncture][1];
        interp_z[ipuncture] = m_puncture_coords[ipuncture][2];
    }

    // setup query
    InterpolationQuery query(m_num_punctures);
    query.setCoords(0, interp_x.data())
        .setCoords(1, interp_y.data())
        .setCoords(2, interp_z.data())
        .addComp(c_shift1, interp_shift1.data())
        .addComp(c_shift2, interp_shift2.data())
        .addComp(c_shift3, interp_shift3.data());

    // engage!
    m_interpolator->interp(query);

    // put the shift values into the output array
    for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
    {
        m_puncture_shift[ipuncture] = {interp_shift1[ipuncture],
                                       interp_shift2[ipuncture],
                                       interp_shift3[ipuncture]};
    }
#endif
    for (int ilevel = 0; ilevel < m_gr_amr->finestLevel(); ++ilevel)
    {
        for (PunctureIter punc_iter(s_particle_container, ilevel);
             punc_iter.isValid(); ++punc_iter)
        {
            const Box &punc_box = punc_iter.validbox();
            const auto &state_mf =
                m_gr_amr->getLevel(ilevel).get_new_data(State_Type);

            const FArrayBox &punc_fab = state_mf[punc_iter];

            amrex::ParallelFor(punc_box, [=] AMREX_GPU_DEVICE())
        }
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
