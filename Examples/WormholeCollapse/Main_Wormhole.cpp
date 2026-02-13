#include "DefaultLevelFactory.hpp"
#include "GRAMR.hpp"
#include "GRParmParse.hpp"
#include "SetupFunctions.hpp"
#include "SimulationParameters.hpp"
#include "WormholeLevel.hpp"

int runGRTeclyn(int /*argc*/, char * /*argv*/[])
{
    BL_PROFILE("runGRTeclyn()");

    GRParmParse pp;
    SimulationParameters sim_params(pp);

    if (sim_params.just_check_params)
        return 0;

    GRAMR::set_simulation_parameters(sim_params);

    // Factory to build Wormhole Levels
    DefaultLevelFactory<WormholeLevel> wh_level_bld;

    // Standard GRAMR container (No need for BHAMR/TPAMR)
    GRAMR wh_amr(&wh_level_bld);

    wh_amr.init(0., sim_params.stop_time);

    while (
        (wh_amr.okToContinue() != 0) &&
        (wh_amr.levelSteps(0) < sim_params.max_steps ||
         sim_params.max_steps < 0) &&
        (wh_amr.cumTime() < sim_params.stop_time || sim_params.stop_time < 0.0))
    {
        wh_amr.coarseTimeStep(sim_params.stop_time);
    }

    if (wh_amr.stepOfLastCheckPoint() < wh_amr.levelSteps(0) &&
        sim_params.checkpoint_interval >= 0)
    {
        wh_amr.checkPoint();
    }

    if (wh_amr.stepOfLastPlotFile() < wh_amr.levelSteps(0) &&
        sim_params.plot_interval >= 0)
    {
        wh_amr.writePlotFile();
    }

    return 0;
}

int main(int argc, char *argv[])
{
    mainSetup(argc, argv);
    int status = runGRTeclyn(argc, argv);
    mainFinalize();
    return status;
}