#include "DefaultLevelFactory.hpp"
#include "GRAMR.hpp"
#include "GRParmParse.hpp"
#include "RadialRecipeLevel.hpp"
#include "RLObservationCollector.hpp"
#include "RadialRecipeConstraintNorms.hpp"
#include "RLActionApplier.hpp"
#include "RLBridge.hpp"
#include "RLMatterPumpParams.hpp"
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

#ifdef USE_RL
    if (sim_params.rl_enabled)
    {
        get_rl_bridge().configure(sim_params.rl_zmq_port,
                                  sim_params.rl_zmq_timeout_ms,
                                  sim_params.rl_zmq_timeout_ms);
        get_rl_bridge().ensure_bound();
        get_rl_bridge().reset_terminate();
    }
#endif

    while (
        (recipe_amr.okToContinue() != 0) &&
        (recipe_amr.levelSteps(0) < sim_params.max_steps ||
         sim_params.max_steps < 0) &&
        (recipe_amr.cumTime() < sim_params.stop_time || sim_params.stop_time < 0.0))
    {
        recipe_amr.coarseTimeStep(sim_params.stop_time);

        if (sim_params.rl_enabled)
        {
            auto &level0 =
                dynamic_cast<RadialRecipeLevel &>(recipe_amr.getLevel(0));
            const auto norms = compute_radial_recipe_constraint_norms(level0);
            RLRuntime::publish_cached_L2_Ham(norms.L2_Ham);

            if (sim_params.rl_coarse_step_interval > 0 &&
                recipe_amr.levelSteps(0) % sim_params.rl_coarse_step_interval ==
                    0)
            {
                const auto obs =
                    collect_rl_observations(recipe_amr, norms.L2_Ham,
                                            norms.L2_Mom);
                std::vector<double> obs_vec(obs.begin(), obs.end());
                auto actions = get_rl_bridge().exchange(
                    obs_vec, 5, sim_params.rl_zmq_timeout_ms);
                std::array<double, 5> action_arr{};
                for (int i = 0; i < 5; ++i)
                {
                    action_arr[static_cast<std::size_t>(i)] = actions[i];
                }
                apply_rl_actions(sim_params.rl_pump_amplitude,
                                 sim_params.rl_pump_frequency,
                                 sim_params.rl_pump_phase,
                                 sim_params.rl_pump_max_amplitude,
                                 sim_params.ccz4_params, action_arr);
            }

#ifdef USE_RL
            if (get_rl_bridge().terminate_requested())
            {
                break;
            }
#endif
        }
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
