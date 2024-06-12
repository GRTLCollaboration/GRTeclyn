/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// Our includes
#include "DefaultLevelFactory.hpp"
#include "GRParmParse.hpp"
#include "MultiLevelTask.hpp"
#include "SetupFunctions.hpp"
#include "SimulationParameters.hpp"

// Problem specific includes:
#include "KleinGordon.hpp" // TPAMR code conditional compiled on USE_TWOPUNCTURES
#include "KleinGordonLevel.hpp"

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

    DefaultLevelFactory<KleinGordonLevel> KleinGordon_level_bld;

    KleinGordon amr(&KleinGordon_level_bld);

    amr.init(0., sim_params.stop_time);

    while ((amr.okToContinue() != 0) &&
           (amr.levelSteps(0) < sim_params.max_steps ||
            sim_params.max_steps < 0) &&
           (amr.cumTime() < sim_params.stop_time || sim_params.stop_time < 0.0))
    {
        amr.coarseTimeStep(sim_params.stop_time);
    }

    // Write final checkpoint and plotfile
    if (amr.stepOfLastCheckPoint() < amr.levelSteps(0) &&
        sim_params.checkpoint_interval >= 0)
    {
        amr.checkPoint();
    }

    if (amr.stepOfLastPlotFile() < amr.levelSteps(0) &&
        sim_params.plot_interval >= 0)
    {
        amr.writePlotFile();
    }

    return 0;
}

int main(int argc, char *argv[])
{
    mainSetup(argc, argv);

    int status = runGRTeclyn(argc, argv);

    if (status == 0)
    {
        amrex::Print() << "GRChombo finished." << std::endl;
    }
    else
    {
        amrex::Print() << "GRChombo failed with return code " << status
                       << std::endl;
    }

    mainFinalize();
    return status;
}
