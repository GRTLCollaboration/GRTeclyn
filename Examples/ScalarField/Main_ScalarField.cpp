/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#include "DefaultLevelBld.hpp"
#include "GRAmr.hpp"
#include "GRParmParse.hpp"
#include "SetupFunctions.hpp"
#include "SimulationParameters.hpp"

#include "ScalarFieldAmr.hpp"
#include "ScalarFieldLevel.hpp"

int runGRTeclyn()
{
    BL_PROFILE("runGRTeclyn()");

    GRParmParse pp; // NOLINT(readability-identifier-length)

    if (just_check_params())
    {
        return 0;
    }

    DefaultLevelBld<ScalarFieldLevel> scalar_field_level_bld;
    ScalarFieldAmr gr_amr(&scalar_field_level_bld);

    amrex::Real stop_time{};
    pp.get("evolution.stop_time", stop_time);
    int max_steps{};
    pp.get("evolution.max_steps", max_steps);

    gr_amr.init(0.0, stop_time);

    // Engage! Run the evolution
    while ((gr_amr.okToContinue() != 0) &&
           (gr_amr.levelSteps(0) < max_steps || max_steps < 0) &&
           (gr_amr.cumTime() < stop_time || stop_time < 0.0))
    {
        gr_amr.coarseTimeStep(stop_time);
    }

    int check_int{};
    pp.get("amr.check_int", check_int);
    int plot_int{};
    pp.get("amr.plot_int", plot_int);

    if (gr_amr.stepOfLastCheckPoint() < gr_amr.levelSteps(0) && check_int >= 0)
    {
        gr_amr.checkPoint();
    }

    if (gr_amr.stepOfLastPlotFile() < gr_amr.levelSteps(0) && plot_int >= 0)
    {
        gr_amr.writePlotFile();
    }

    return 0;
}

int main(int argc, char *argv[])
{
    mainSetup(argc, argv);

    const int status = runGRTeclyn();

    print_job_end_message(status);

    mainFinalize();
    return status;
}
