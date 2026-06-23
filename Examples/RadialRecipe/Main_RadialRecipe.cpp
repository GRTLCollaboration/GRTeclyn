#include "DefaultLevelFactory.hpp"
#include "GRAMR.hpp"
#include "GRParmParse.hpp"
#include "RadialRecipeLevel.hpp"
#include "SetupFunctions.hpp"
#include "SimulationParameters.hpp"

#ifdef USE_RL
#include "RLObservationCollector.hpp"
#include "RadialRecipeConstraintNorms.hpp"
#include "RLActionApplier.hpp"
#include "RLBridge.hpp"
#include "RLMatterPumpParams.hpp"
#endif

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
        RLRuntime::seed_lumps(sim_params.rl_num_lumps,
                              sim_params.rl_lump_seed_x.data(),
                              sim_params.rl_lump_seed_y.data(),
                              sim_params.rl_lump_seed_z.data());
    }
#endif

    while (
        (recipe_amr.okToContinue() != 0) &&
        (recipe_amr.levelSteps(0) < sim_params.max_steps ||
         sim_params.max_steps < 0) &&
        (recipe_amr.cumTime() < sim_params.stop_time || sim_params.stop_time < 0.0))
    {
        recipe_amr.coarseTimeStep(sim_params.stop_time);

#ifdef USE_RL
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
                const bool complex_field =
                    (sim_params.recipe_matter_model ==
                     "grtresna_complex_scalar");
                const double ball_radius = 3.0 * sim_params.rl_pump_width;
                const int num_lumps      = sim_params.rl_num_lumps;
                const auto obs           = collect_rl_observations(
                    recipe_amr, norms.L2_Ham, norms.L2_Mom, num_lumps,
                    complex_field, ball_radius,
                    sim_params.recipe_params.grid_center);
                const int action_dim = 3 * num_lumps + 2;
                auto actions         = get_rl_bridge().exchange(
                    obs, action_dim, sim_params.rl_zmq_timeout_ms);
                apply_rl_actions(num_lumps, sim_params.rl_pump_amplitude,
                                 sim_params.rl_pump_frequency,
                                 sim_params.rl_pump_phase,
                                 sim_params.rl_pump_max_amplitude,
                                 sim_params.ccz4_params, actions);
            }

            if (get_rl_bridge().terminate_requested())
            {
                break;
            }
        }
#endif
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
