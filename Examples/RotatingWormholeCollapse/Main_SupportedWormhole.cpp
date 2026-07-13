#include "BHAMR.hpp"
#include "DefaultLevelFactory.hpp"
#include "GRAMR.hpp"
#include "GRParmParse.hpp"
#include "RLMatterPumpParams.hpp"
#include "SetupFunctions.hpp"
#include "SimulationParameters.hpp"
#include "SupportedWormholeLevel.hpp"
#include "TrajectoryEvaluator.hpp"

#include <array>

//! Pump ramp factor: 1 before t_start, linear from 1→floor over
//! [t_start, t_end], floor after.  t_start < 0 => always 1 (no ramp).
inline double pump_ramp_factor(double t, double t_start, double t_end,
                               double floor)
{
    if (t_start < 0.0 || t < t_start)
        return 1.0;
    if (t_end <= t_start || t >= t_end)
        return floor;
    const double s = (t - t_start) / (t_end - t_start);
    return 1.0 + s * (floor - 1.0);
}

int runGRTeclyn(int /*argc*/, char * /*argv*/[])
{
    BL_PROFILE("runGRTeclyn()");

    GRParmParse pp;
    SimulationParameters sim_params(pp);

    if (sim_params.just_check_params)
        return 0;

    GRAMR::set_simulation_parameters(sim_params);

    DefaultLevelFactory<SupportedWormholeLevel> wh_level_bld;

    // BHAMR (a GRAMR child) so that BHAMR::init() sets up the Weyl4
    // ParticleInterpolator used for in-code Psi4 spherical-harmonic extraction
    // (see SupportedWormholeLevel::specificPostTimeStep).  num_punctures is a
    // no-op here; puncture tracking is disabled.
    BHAMR<SupportedWormholeLevel::num_punctures> wh_amr(&wh_level_bld);

    wh_amr.init(0., sim_params.stop_time);

    // --- Trajectory-guided matter pump (Rung 2 active support) ---------------
    // When trajectory_mode == 1, pump site centres are driven by parametric
    // orbit equations evaluated on the CPU each coarse step.  No ZMQ/RL needed.
    if (sim_params.trajectory_mode == 1)
    {
        const auto &tp = sim_params.trajectory_params;
        const int n_traj = tp.num_lumps;
        std::array<double, RL_MAX_LUMPS> tx{}, ty{}, tz{}, tamp{};
        TrajectoryEvaluator::evaluate(
            0.0, tp, tx.data(), ty.data(), tz.data(), tamp.data());
        RLRuntime::seed_lumps(n_traj, tx.data(), ty.data(), tz.data());
        for (int k = 0; k < n_traj; ++k)
        {
            sim_params.rl_pump_amplitude[k] = tamp[k];
            sim_params.rl_pump_frequency[k] =
                sim_params.trajectory_pump_frequency;
        }
        sim_params.rl_pump_width = tp.well_width;
    }

    // Base exotic-support coupling; the support ramp (Phase 6 collapse trigger)
    // scales this each coarse step.  Captured before the loop so the ramp is
    // relative to the configured value, not the previous step's ramped value.
    const double base_support_strength =
        sim_params.wormhole_params.support_strength;

    while (
        (wh_amr.okToContinue() != 0) &&
        (wh_amr.levelSteps(0) < sim_params.max_steps ||
         sim_params.max_steps < 0) &&
        (wh_amr.cumTime() < sim_params.stop_time || sim_params.stop_time < 0.0))
    {
        wh_amr.coarseTimeStep(sim_params.stop_time);

        // --- Support-strength ramp (Phase 6 collapse-on-command) -----------
        // Independent of the pump: ramp the exotic stress-energy coupling from
        // base -> base*floor over [t_start, t_end].  simParams() is a live
        // reference (GRAMR stores &sim_params), so the matter class picks up
        // the new value on the next RHS eval.
        {
            const double t = wh_amr.cumTime();
            const double sramp = pump_ramp_factor(
                t, sim_params.support_ramp_t_start,
                sim_params.support_ramp_t_end, sim_params.support_ramp_floor);
            sim_params.wormhole_params.support_strength =
                base_support_strength * sramp;
        }

        // --- Trajectory update (CPU, every coarse step) --------------------
        if (sim_params.trajectory_mode == 1)
        {
            const double t   = wh_amr.cumTime();
            const auto &tp   = sim_params.trajectory_params;
            const int n_traj = tp.num_lumps;
            std::array<double, RL_MAX_LUMPS> tx{}, ty{}, tz{}, tamp{};
            TrajectoryEvaluator::evaluate(
                t, tp, tx.data(), ty.data(), tz.data(), tamp.data());
            const double ramp = pump_ramp_factor(
                t, sim_params.pump_ramp_t_start, sim_params.pump_ramp_t_end,
                sim_params.pump_ramp_floor);
            for (int k = 0; k < n_traj; ++k)
            {
                RLRuntime::g_lump_state[k].x = tx[k];
                RLRuntime::g_lump_state[k].y = ty[k];
                RLRuntime::g_lump_state[k].z = tz[k];
                sim_params.rl_pump_amplitude[k] = tamp[k] * ramp;
            }

            // Publish L2_Ham for the governor.  The constraint norms are
            // computed in specificPostTimeStep (Level 0) and written to
            // constraint_norms.dat; read the latest value from there.
            // (The wormhole example doesn't have a constraint-norms helper
            // like RadialRecipe, so we read the file.  This is cheap — one
            // file read per coarse step.)
            // For now, use the cached value (set by specificPostTimeStep if
            // it publishes; otherwise the governor defaults to 1.0).
        }
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