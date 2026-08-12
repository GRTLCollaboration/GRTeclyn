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

#include "ParticleInterpolator.hpp"

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

    if (pp.contains("check_params"))
    {
        return 0;
    }

    DefaultLevelFactory<BinaryBHLevel> bh_level_bld;

#ifdef USE_TWOPUNCTURES
    TPAMR bh_amr;
    bh_amr.set_two_punctures_parameters(sim_params.tp_params);
    // Run TwoPunctures solver
    bh_amr.m_two_punctures.Run();
#else
    BHAMR<BinaryBHLevel::num_punctures> bh_amr(&bh_level_bld);
#endif

    double stop_time{};
    pp.get("amr.stop_time", stop_time);
    int max_steps{};
    pp.get("amr.max_steps", max_steps);

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
