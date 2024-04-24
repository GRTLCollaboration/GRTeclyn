/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// Our includes
#include "DefaultLevelFactory.hpp"
#include "GRParmParse.hpp"
#include "MultiLevelTask.hpp"
#include "SetupFunctions.hpp"
#include "SimulationParameters.hpp"
// TPAMR.hpp includes BHAMR.hpp
#include "TPAMR.hpp" // TPAMR code conditional compiled on USE_TWOPUNCTURES

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

#ifdef USE_TWOPUNCTURES
    TPAMR bh_amr;
    bh_amr.set_two_punctures_parameters(sim_params.tp_params);
    // Run TwoPunctures solver
    bh_amr.m_two_punctures.Run();
#else
    BHAMR bh_amr(&bh_level_bld);
#endif

    bh_amr.init(0., sim_params.stop_time);

    auto start_time = std::chrono::steady_clock::now();

    while (
        (bh_amr.okToContinue() != 0) &&
        (bh_amr.levelSteps(0) < sim_params.max_steps ||
         sim_params.max_steps < 0) &&
        (bh_amr.cumTime() < sim_params.stop_time || sim_params.stop_time < 0.0))
    {
        bh_amr.coarseTimeStep(sim_params.stop_time);
    }

    auto end_time = std::chrono::steady_clock::now();
    auto elapsed  = std::chrono::duration_cast<std::chrono::duration<double>>(
        end_time - start_time);

    amrex::Print().SetPrecision(16)
        << "Total simulation time = " << elapsed.count() << " secs\n";

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
        amrex::Print() << "GRTeclyn finished." << std::endl;
    }
    else
    {
        amrex::Print() << "GRTeclyn failed with return code " << status
                       << std::endl;
    }

    mainFinalize();
    return status;
}
