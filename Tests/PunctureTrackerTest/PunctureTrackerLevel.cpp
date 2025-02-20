/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// GRTeclyn and test headers
#include "PunctureTrackerLevel.hpp"
#include "PunctureTagger.hpp"
#include "PunctureTracker.hpp"

// doctest header
#include "doctest.h"

BHAMR<PunctureTrackerLevel::num_punctures> *
PunctureTrackerLevel::get_bhamr_ptr()
{
    return dynamic_cast<BHAMR<num_punctures> *>(get_gramr_ptr());
}

PunctureTracker<PunctureTrackerLevel::num_punctures> &
PunctureTrackerLevel::get_puncture_tracker()
{
    return get_bhamr_ptr()->get_puncture_tracker();
}

void PunctureTrackerLevel::variableSetUp() { stateVariableSetUp(); }

void PunctureTrackerLevel::initData()
{
    amrex::MultiFab &state = get_new_data(State_Type);
    const auto &arrs       = state.arrays();
    amrex::ParallelFor(state, state.nGrowVect(),
                       // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       {
                           const auto &array = arrs[box_no];
                           for (int icomp = 0; icomp < array.nComp(); ++icomp)
                           {
                               array(i, j, k, icomp) = 0.0;
                           }
                           array(i, j, k, c_shift2) = shift_y_val;
                       });
    amrex::Gpu::streamSynchronize();

    if (Level() == 0)
    {
        // need to set the puncture coordinates as we use it for the puncture
        // tagging
        get_puncture_tracker().set_puncture_coords(
            simParams().puncture_tracking_initial_coords);
        // can't call start_from_initial_punctures() because we need the full
        // AMR grid first
    }
}

// Calculate RHS during RK4 substeps
// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
void PunctureTrackerLevel::specificEvalRHS(amrex::MultiFab &a_soln,
                                           amrex::MultiFab &a_rhs,
                                           const double /*a_time*/)
{
    // We don't need any evolution in this test.
    a_rhs.setVal(0.0);
}

void PunctureTrackerLevel::tag_cells(amrex::TagBoxArray &a_tag_box_array,
                                     amrex::Real a_regrid_threshold)
{
    amrex::MultiFab &state_new = get_new_data(State_Type);

    const auto &tag_arrs       = a_tag_box_array.arrays();
    const auto &state_new_arrs = state_new.const_arrays();

    std::array<amrex::Real, num_puncture_coords> puncture_coords{};

    puncture_coords = get_puncture_tracker().get_puncture_coords();

    PunctureTagger<num_punctures> puncture_tagger(
        Geom().CellSize(0), Level(), get_gramr_ptr()->maxLevel(),
        puncture_coords,
        {simParams().fake_bh1_mass, simParams().fake_bh2_mass});

    amrex::ParallelFor(state_new, amrex::IntVect(0),
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       { puncture_tagger(i, j, k, tag_arrs[box_no]); });
    amrex::Gpu::streamSynchronize();
}

void PunctureTrackerLevel::specific_post_init()
{
    get_puncture_tracker().start_from_initial_punctures();
}

void PunctureTrackerLevel::specific_post_restart()
{

    std::string restart_checkpoint{};
    GRParmParse pp("amr");
    pp.get("restart", restart_checkpoint);
    get_puncture_tracker().restart(restart_checkpoint);
}

void PunctureTrackerLevel::specific_post_checkpoint(
    const std::string &a_chk_dir, std::ostream & /*a_os*/)
{
    get_puncture_tracker().checkpoint(a_chk_dir);
}

void PunctureTrackerLevel::specificPostTimeStep()
{
    if (Level() == simParams().puncture_tracking_level)
    {
        bool write_punctures = false;
        amrex::Real cur_time = get_state_data(State_Type).curTime();
        amrex::Real dt       = get_gramr_ptr()->dtLevel(Level());
        get_puncture_tracker().track(cur_time, dt, write_punctures);

        auto correct_puncture_coords =
            simParams().puncture_tracking_initial_coords;
        for (int ipuncture = 0; ipuncture < num_punctures; ++ipuncture)
        {
            correct_puncture_coords[ipuncture * AMREX_SPACEDIM + 1] -=
                cur_time * shift_y_val;
        }
        auto puncture_coords = get_puncture_tracker().get_puncture_coords();
        constexpr int cout_precision = 16;
        for (int icoord = 0; icoord < num_puncture_coords; ++icoord)
        {
            INFO("puncture_coords["
                 << icoord << "] = " << std::setprecision(cout_precision)
                 << puncture_coords[icoord]);
            INFO("correct_puncture_coords["
                 << icoord << "] = " << std::setprecision(cout_precision)
                 << correct_puncture_coords[icoord]);
            CHECK(puncture_coords[icoord] ==
                  doctest::Approx(correct_puncture_coords[icoord])
                      .epsilon(1e-12));
        }
    }
}
