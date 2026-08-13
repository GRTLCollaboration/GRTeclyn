/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#include "DefaultLevelFactory.hpp"
#include "GRAMR.hpp"
#include "GRParmParse.hpp"
#include "SetupFunctions.hpp"
#include "SimulationParameters.hpp"

#include "ScalarFieldAMR.hpp"
#include "ScalarFieldLevel.hpp"

// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
int runGRTeclyn(int /*argc*/, char * /*argv*/[])
{
    BL_PROFILE("runGRTeclyn()");

    GRParmParse pp; // NOLINT(readability-identifier-length)
    SimulationParameters sim_params(pp);

    if (sim_params.just_check_params)
    {
        return 0;
    }

    GRAMR::set_simulation_parameters(sim_params);
    DefaultLevelFactory<ScalarFieldLevel> scalar_field_level_bld;
    ScalarFieldAMR gr_amr(&scalar_field_level_bld);

    gr_amr.init(0.0, sim_params.stop_time);

    while (
        (gr_amr.okToContinue() != 0) &&
        (gr_amr.levelSteps(0) < sim_params.max_steps ||
         sim_params.max_steps < 0) &&
        (gr_amr.cumTime() < sim_params.stop_time || sim_params.stop_time < 0.0))
    {
        gr_amr.coarseTimeStep(sim_params.stop_time);
    }

    if (gr_amr.stepOfLastCheckPoint() < gr_amr.levelSteps(0) &&
        sim_params.checkpoint_interval >= 0)
    {
        gr_amr.checkPoint();
    }

    if (gr_amr.stepOfLastPlotFile() < gr_amr.levelSteps(0) &&
        sim_params.plot_interval >= 0)
    {
        gr_amr.writePlotFile();
    }

    return 0;
}

int main(int argc, char *argv[])
{
    mainSetup(argc, argv);

    const int status = runGRTeclyn(argc, argv);

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
