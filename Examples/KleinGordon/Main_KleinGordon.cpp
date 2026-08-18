/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// Our includes
#include "DefaultLevelFactory.hpp"
#include "GRParmParse.hpp"
#include "SetupFunctions.hpp"
#include "SimulationParameters.hpp"

// Problem specific includes:
#include "KleinGordonLevel.hpp"

// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
int runGRTeclyn(int argc, char *argv[])
{
    BL_PROFILE("runGRTeclyn()");

    // Load the parameter file and construct the SimulationParameter class
    // To add more parameters edit the SimulationParameters file.
    GRParmParse pp; // NOLINT(readability-identifier-length)

    if (just_check_params(argc, argv))
    {
        return 0;
    }

    std::string model{};
    pp.get("klein_gordon.model", model);

    double stop_time{};
    pp.get("evolution.stop_time", stop_time);
    int max_steps{};
    pp.get("evolution.max_steps", max_steps);

    amrex::Print() << "Now running " << model << " simulation"
                   << "\n";

    DefaultLevelFactory<KleinGordonLevel> KleinGordon_level_bld;

    GRAMR amr(&KleinGordon_level_bld);

    amr.init(0., stop_time);

    while ((amr.okToContinue() != 0) &&
           (amr.levelSteps(0) < max_steps || max_steps < 0) &&
           (amr.cumTime() < stop_time || stop_time < 0.0))
    {
        amr.coarseTimeStep(stop_time);
    }

    int check_int{}; // Steps between checkpoint file outputs
    pp.get("amr.check_int", check_int);

    int plot_int{}; // Steps between plot file outputs
    pp.get("amr.plot_int", plot_int);

    // Write final checkpoint and plotfile
    if (amr.stepOfLastCheckPoint() < amr.levelSteps(0) && check_int >= 0)
    {
        amr.checkPoint();
    }

    if (amr.stepOfLastPlotFile() < amr.levelSteps(0) && plot_int >= 0)
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
        amrex::Print() << "GRTeclyn finished.\n";
    }
    else
    {
        amrex::Print() << "GRTeclyn failed with return code " << status << "\n";
    }

    mainFinalize();
    return status;
}
