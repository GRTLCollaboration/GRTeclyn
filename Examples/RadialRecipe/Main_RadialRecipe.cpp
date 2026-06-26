#include "DefaultLevelFactory.hpp"
#include "GRAMR.hpp"
#include "GRParmParse.hpp"
#include "RadialRecipeLevel.hpp"
#include "SetupFunctions.hpp"
#include "SimulationParameters.hpp"

#include "TrajectoryEvaluator.hpp"
#include "TrajectoryParams.hpp"
#include "RadialRecipeConstraintNorms.hpp"
#include "RLMatterPumpParams.hpp"

#ifdef USE_RL
#include "RLObservationCollector.hpp"
#include "RLActionApplier.hpp"
#include "RLBridge.hpp"
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

    // --- Trajectory-guided geometry survey (no RL / no ZMQ) -----------------
    // When trajectory_mode == 1, the pump site centres are driven by parametric
    // equations evaluated on the CPU each coarse step.  The pump amplitudes come
    // from the trajectory params (not from an RL agent), so rl_enabled is not
    // required.  We seed the lump tracker so the pump struct is well-formed even
    // if no RL bridge is active.
    if (sim_params.trajectory_mode == 1)
    {
        const auto &tp = sim_params.trajectory_params;
        const int n_traj = tp.num_lumps;
        std::array<double, GRTRESNA_MAX_INDEPENDENT_SCALARS> tx{}, ty{}, tz{};
        std::array<double, GRTRESNA_MAX_INDEPENDENT_SCALARS> tamp{};
        TrajectoryEvaluator::evaluate(
            0.0, tp, tx.data(), ty.data(), tz.data(), tamp.data());
        // Seed lump state with t=0 trajectory positions.
        RLRuntime::seed_lumps(n_traj, tx.data(), ty.data(), tz.data());
        // Set pump amplitudes and frequency from trajectory params.
        for (int k = 0; k < n_traj; ++k)
        {
            sim_params.rl_pump_amplitude[k] = tamp[k];
            sim_params.rl_pump_frequency[k] =
                sim_params.trajectory_pump_frequency;
        }
        sim_params.rl_num_lumps  = n_traj;
        sim_params.rl_pump_width = tp.well_width;
    }

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

        // --- Trajectory update (CPU, every coarse step) --------------------
        // Evaluate parametric equations at current simulation time.
        // Updates lump centres and pump amplitudes in sim_params so the next
        // RHS evaluation uses the new positions.  No GPU trig — only fast
        // FMA distance calculations in the kernel.
        if (sim_params.trajectory_mode == 1)
        {
            const double t   = recipe_amr.cumTime();
            const auto &tp   = sim_params.trajectory_params;
            const int n_traj = tp.num_lumps;
            std::array<double, GRTRESNA_MAX_INDEPENDENT_SCALARS> tx{}, ty{},
                tz{};
            std::array<double, GRTRESNA_MAX_INDEPENDENT_SCALARS> tamp{};
            TrajectoryEvaluator::evaluate(
                t, tp, tx.data(), ty.data(), tz.data(), tamp.data());
            // Write directly into global lump state (pump builder reads these).
            for (int k = 0; k < n_traj; ++k)
            {
                RLRuntime::g_lump_state[k].x = tx[k];
                RLRuntime::g_lump_state[k].y = ty[k];
                RLRuntime::g_lump_state[k].z = tz[k];
                sim_params.rl_pump_amplitude[k] = tamp[k];
            }

            // Publish L2 Ham for the governor (even without RL).
            auto &level0 =
                dynamic_cast<RadialRecipeLevel &>(recipe_amr.getLevel(0));
            const auto norms =
                compute_radial_recipe_constraint_norms(level0);
            RLRuntime::publish_cached_L2_Ham(norms.L2_Ham);
        }

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
                         "grtresna_complex_scalar" ||
                     sim_params.recipe_matter_model ==
                         "grtresna_bicomplex_scalar");
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
