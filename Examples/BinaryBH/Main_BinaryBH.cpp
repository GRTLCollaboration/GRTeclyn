/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// Our includes
#include "BHAmr.hpp"
#include "DefaultLevelBld.hpp"
#include "GRParmParse.hpp"
#include "MultiLevelTask.hpp"
#include "SetupFunctions.hpp"
#include "SimulationParameters.hpp"

// Problem specific includes:
#include "BinaryBHLevel.hpp"

// System includes
#include <chrono>
#include <iostream>

int runGRTeclyn()
{
    BL_PROFILE("runGRTeclyn()");

    // Load the parameter file and construct the SimulationParameter class
    // To add more parameters edit the SimulationParameters file.
    GRParmParse pp; // NOLINT(readability-identifier-length)

    if (just_check_params())
    {
        return 0;
    }

    DefaultLevelBld<BinaryBHLevel> bh_level_bld;

    BHAmr<BinaryBHLevel::num_punctures> bh_amr(&bh_level_bld);

    double stop_time{};
    pp.get("evolution.stop_time", stop_time);
    int max_steps{};
    pp.get("evolution.max_steps", max_steps);

    bh_amr.init(0., stop_time);

    while ((bh_amr.okToContinue() != 0) &&
           (bh_amr.levelSteps(0) < max_steps || max_steps < 0) &&
           (bh_amr.cumTime() < stop_time || stop_time < 0.0))
    {
        bh_amr.coarseTimeStep(stop_time);
    }

    int check_int{}; // Steps between checkpoint file outputs
    pp.get("amr.check_int", check_int);

    int plot_int{}; // Steps between plot file outputs
    pp.get("amr.plot_int", plot_int);

    // Write final checkpoint and plotfile
    if (bh_amr.stepOfLastCheckPoint() < bh_amr.levelSteps(0) && check_int >= 0)
    {
        bh_amr.checkPoint();
    }

    if (bh_amr.stepOfLastPlotFile() < bh_amr.levelSteps(0) && plot_int >= 0)
    {
        bh_amr.writePlotFile();
    }

    return 0;
}

int main(int argc, char *argv[])
{
    mainSetup(argc, argv);

    int status = runGRTeclyn();

    print_job_end_message(status);

    mainFinalize();
    return status;
}
