/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// Our includes
#include "BHAMR.hpp"
#include "DefaultLevelFactory.hpp"
#include "GRParmParse.hpp"
#include "LIKWIDMarkers.hpp"
#include "MultiLevelTask.hpp"
#include "SetupFunctions.hpp"
#include "SimulationParameters.hpp"

// Problem specific includes:
#include "BinaryBHLevel.hpp"

// System includes
#include <chrono>
#include <iostream>

// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
int runGRTeclyn(int /*argc*/, char * /*argv*/[])
{
    BL_PROFILE("runGRTeclyn()");

    // Load the parameter file and construct the SimulationParameter class
    // To add more parameters edit the SimulationParameters file.
    GRParmParse pp; // NOLINT(readability-identifier-length)
    SimulationParameters sim_params(pp);

    if (sim_params.just_check_params)
    {
        return 0;
    }

    GRAMR::set_simulation_parameters(sim_params);
    DefaultLevelFactory<BinaryBHLevel> bh_level_bld;

    BHAMR<BinaryBHLevel::num_punctures> bh_amr(&bh_level_bld);

    bh_amr.init(0., sim_params.stop_time);

    LIKWID_MARKER_INIT;
    while (
        (bh_amr.okToContinue() != 0) &&
        (bh_amr.levelSteps(0) < sim_params.max_steps ||
         sim_params.max_steps < 0) &&
        (bh_amr.cumTime() < sim_params.stop_time || sim_params.stop_time < 0.0))
    {
        bh_amr.coarseTimeStep(sim_params.stop_time);
    }
    LIKWID_MARKER_CLOSE;

    // Write final checkpoint and plotfile
    if (bh_amr.stepOfLastCheckPoint() < bh_amr.levelSteps(0) &&
        sim_params.checkpoint_interval >= 0)
    {
        bh_amr.checkPoint();
    }

    if (bh_amr.stepOfLastPlotFile() < bh_amr.levelSteps(0) &&
        sim_params.plot_interval >= 0)
    {
        bh_amr.writePlotFile();
    }

    return 0;
}

int main(int argc, char *argv[])
{
    mainSetup(argc, argv);

    int status = runGRTeclyn(argc, argv);

    if (status == 0)
    {
        amrex::Print() << "GRTeclyn finished."
                       << "\n";
    }
    else
    {
        amrex::Print() << "GRTeclyn failed with return code " << status << "\n";
    }

    mainFinalize();
    return status;
}
