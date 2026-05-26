#include "DefaultLevelFactory.hpp"
#include "GRAMR.hpp"
#include "GRParmParse.hpp"
#include "RadialRecipeLevel.hpp"
#include "SetupFunctions.hpp"
#include "SimulationParameters.hpp"

int runGRTeclyn(int /*argc*/, char * /*argv*/[])
{
    BL_PROFILE("runGRTeclyn()");

    GRParmParse pp;
    SimulationParameters sim_params(pp);

    if (sim_params.just_check_params)
        return 0;

    GRAMR::set_simulation_parameters(sim_params);

    DefaultLevelFactory<RadialRecipeLevel> recipe_level_bld;

    GRAMR recipe_amr(&recipe_level_bld);

    recipe_amr.init(0., sim_params.stop_time);

    while (
        (recipe_amr.okToContinue() != 0) &&
        (recipe_amr.levelSteps(0) < sim_params.max_steps ||
         sim_params.max_steps < 0) &&
        (recipe_amr.cumTime() < sim_params.stop_time || sim_params.stop_time < 0.0))
    {
        recipe_amr.coarseTimeStep(sim_params.stop_time);
    }

    if (recipe_amr.stepOfLastCheckPoint() < recipe_amr.levelSteps(0) &&
        sim_params.checkpoint_interval >= 0)
    {
        recipe_amr.checkPoint();
    }

    if (recipe_amr.stepOfLastPlotFile() < recipe_amr.levelSteps(0) &&
        sim_params.plot_interval >= 0)
    {
        recipe_amr.writePlotFile();
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
