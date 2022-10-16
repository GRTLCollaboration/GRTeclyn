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
// TPAMR.hpp includes BHAMR.hpp
#include "TPAMR.hpp" // TPAMR code conditional compiled on USE_TWOPUNCTURES

// Problem specific includes:
#include "BinaryBHLevel.hpp"

// System includes
#include <chrono>
#include <iostream>

DefaultLevelFactory<BinaryBHLevel> bh_level_bld;

int runGRAMReX(int /*argc*/, char */*argv*/[])
{
    BL_PROFILE("runGRAMReX()");

    // Load the parameter file and construct the SimulationParameter class
    // To add more parameters edit the SimulationParameters file.
    GRParmParse pp;
    SimulationParameters sim_params(pp);

    if (sim_params.just_check_params)
        return 0;

    GRAMR::set_simulation_parameters(sim_params);

#ifdef USE_TWOPUNCTURES
    TPAMR bh_amr;
    bh_amr.set_two_punctures_parameters(sim_params.tp_params);
    // Run TwoPunctures solver
    bh_amr.m_two_punctures.Run();
#else
    BHAMR bh_amr(&bh_level_bld);
#endif

    bh_amr.init(0., sim_params.stop_time);

    while (
        bh_amr.okToContinue() &&
        (bh_amr.levelSteps(0) < sim_params.max_steps ||
         sim_params.max_steps < 0) &&
        (bh_amr.cumTime() < sim_params.stop_time || sim_params.stop_time < 0.0))
    {
        bh_amr.coarseTimeStep(sim_params.stop_time);
    }

    // Write final checkpoint and plotfile
    if (bh_amr.stepOfLastCheckPoint() < bh_amr.levelSteps(0))
    {
        bh_amr.checkPoint();
    }

    if (bh_amr.stepOfLastPlotFile() < bh_amr.levelSteps(0))
    {
        bh_amr.writePlotFile();
    }

    return 0;
}

int main(int argc, char *argv[])
{
    mainSetup(argc, argv);

    int status = runGRAMReX(argc, argv);

    if (status == 0)
        amrex::Print() << "GRChombo finished." << std::endl;
    else
        amrex::Print() << "GRChombo failed with return code " << status << std::endl;

    mainFinalize();
    return status;
}
