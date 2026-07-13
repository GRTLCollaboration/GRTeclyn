#ifndef SIMULATIONPARAMETERS_HPP
#define SIMULATIONPARAMETERS_HPP

#include "ExternalGridInitialData.hpp"
#include "GRParmParse.hpp"
#include "RLLumpState.hpp" // RL_MAX_LUMPS
#include "SimulationParametersBase.hpp"
#include "SpongeZone.hpp" // SpongeZoneParams
#include "SupportedWormholeInitialData.hpp"
#include "TrajectoryParams.hpp"

#include <array>
#include <sstream>
#include <string>

class SimulationParameters : public SimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        read_shared_params(pp);
        read_wormhole_params(pp);
        check_params();
    }

    void read_shared_params(GRParmParse &pp)
    {
        pp.load("calculate_constraint_norms", calculate_constraint_norms,
                false);
        pp.load("wormhole_matter_model", wormhole_matter_model,
                std::string("exotic_scalar"));
    }

    void read_wormhole_params(GRParmParse &pp)
    {
        pp.load("wormhole_initial_lapse_type", wormhole_params.initial_lapse_type,
                0);

        pp.load("center", wormhole_params.grid_center, center);

        double b0_single = 1.0;
        std::array<double, AMREX_SPACEDIM> c_single = {0.0, 0.0, 0.0};
        pp.load("wormhole_throat_radius", b0_single, 1.0);
        pp.load("wormhole_centerA", wormhole_params.centerA, c_single);
        wormhole_params.b0 = b0_single;

        pp.load("phantom_mass", wormhole_params.phantom_mass, 0.0);
        // Q-ball self-interaction couplings (see params_t docs).  Must match the
        // GRTresna solve's scalar_lambda / scalar_mu for the loaded .gridinit to
        // stay in equilibrium at t=0.  Default 0 = free field (backward compat).
        pp.load("phantom_lambda", wormhole_params.phantom_lambda, 0.0);
        pp.load("phantom_mu6", wormhole_params.phantom_mu6, 0.0);

        pp.load("wormhole_support_strength", wormhole_params.support_strength,
                1.0);

        pp.load("wormhole_phi_monopole_amplitude",
                wormhole_params.phi_monopole_amplitude, 0.0);
        pp.load("wormhole_phi_perturbation_amplitude",
                wormhole_params.phi_perturbation_amplitude, 0.0);
        pp.load("wormhole_phi_perturbation_width",
                wormhole_params.phi_perturbation_width, 0.0);

        // Spinning complex phantom scalar initial data. Enabled automatically
        // when the matter model is complex_scalar; m and omega set the phase
        // winding and rotation rate.
        wormhole_params.complex_scalar_init =
            (wormhole_matter_model == "complex_scalar") ? 1 : 0;
        pp.load("wormhole_azimuthal_m", wormhole_params.azimuthal_m, 1);
        pp.load("wormhole_rotation_omega", wormhole_params.rotation_omega, 0.0);

        pp.load("recipe_initial_data_file", recipe_initial_data_file,
                std::string(""));
        if (!recipe_initial_data_file.empty())
        {
            external_grid_params.gridinit_file = recipe_initial_data_file;
            external_grid_params.grid_center = center;
        }

        // ---- Trajectory-guided matter pump (Rung 2 active support) ------
        // Opt-in; trajectory_mode=0 (default) is a complete no-op so existing
        // passive runs are bit-identical.  When active, per-lump parametric
        // orbits drive the pump site centres each coarse step (no ZMQ/RL).
        read_trajectory_params(pp);

        // Pump gains / limits (mirror RadialRecipe).  Used only when
        // trajectory_mode == 1.
        pp.load("rl_pump_width", rl_pump_width, 1.5);
        pp.load("rl_pump_max_amplitude", rl_pump_max_amplitude, 0.05);
        pp.load("rl_l2_ham_governor_center", rl_l2_ham_governor_center, 0.035);
        pp.load("rl_l2_ham_governor_width", rl_l2_ham_governor_width, 0.003);
        pp.load("rl_pump_kp", rl_pump_kp, 0.0);
        pp.load("rl_pump_kd", rl_pump_kd, 0.0);
        pp.load("rl_pump_target_profile", rl_pump_target_profile, 0);
        pp.load("rl_pump_target_width", rl_pump_target_width, 0.0);
        pp.load("rl_pump_target_amp", rl_pump_target_amp, 0.0);
        rl_pump_amplitude.fill(0.0);
        rl_pump_frequency.fill(0.0);
        rl_pump_phase.fill(0.0);

        // Pump ramp schedule (Phase 6 collapse trigger).  t_start < 0 => never.
        pp.load("pump_ramp_t_start", pump_ramp_t_start, -1.0);
        pp.load("pump_ramp_t_end", pump_ramp_t_end, 0.0);
        pp.load("pump_ramp_floor", pump_ramp_floor, 0.0);

        // Support-strength ramp schedule (Phase 6 collapse-on-command trigger).
        // Ramps wormhole_support_strength (the exotic stress-energy coupling)
        // from its base value to base*floor over [t_start, t_end] -- the
        // rotating analogue of the static Ellis-Bronnikov S_support cut.
        // t_start < 0 => never (support held constant, backward compatible).
        pp.load("support_ramp_t_start", support_ramp_t_start, -1.0);
        pp.load("support_ramp_t_end", support_ramp_t_end, 0.0);
        pp.load("support_ramp_floor", support_ramp_floor, 0.0);

        // Numerical sponge zone: radially-ramped extra KO dissipation in an
        // outer shell to absorb outgoing waves before the Sommerfeld boundary
        // (clean GW extraction on a larger box).  Disabled by default.
        pp.load("sponge_enabled", sponge_params.enabled, false);
        pp.load("sponge_inner_radius", sponge_params.inner_radius, 24.0);
        pp.load("sponge_outer_radius", sponge_params.outer_radius, 32.0);
        pp.load("sponge_strength", sponge_params.strength, 4.0);
        pp.load("sponge_ramp_power", sponge_params.ramp_power, 4);
        pp.load("sponge_center", sponge_params.center, center);
    }

    void read_trajectory_params(GRParmParse &pp)
    {
        pp.load("trajectory_mode", trajectory_mode, 0);
        {
            int n_traj = 1;
            pp.load("trajectory_num_lumps", n_traj, 1);
            if (n_traj < 1)
                n_traj = 1;
            if (n_traj > RL_MAX_LUMPS)
                n_traj = RL_MAX_LUMPS;
            trajectory_params.num_lumps = n_traj;
        }
        pp.load("trajectory_A_breath", trajectory_params.A_breath, 0.0);
        pp.load("trajectory_omega_breath", trajectory_params.omega_breath, 0.0);
        pp.load("trajectory_z_amp", trajectory_params.z_amp, 0.0);
        pp.load("trajectory_omega_z", trajectory_params.omega_z, 0.0);
        pp.load("trajectory_well_width", trajectory_params.well_width, 1.5);
        pp.load("trajectory_pump_frequency", trajectory_pump_frequency, 0.0);
        for (int k = 0; k < trajectory_params.num_lumps; ++k)
        {
            load_trajectory_lump(pp, k, trajectory_params.lumps[k]);
        }
    }

    void check_params()
    {
        check_parameter("wormhole_initial_lapse_type",
                        wormhole_params.initial_lapse_type,
                        (wormhole_params.initial_lapse_type >= 0) &&
                            (wormhole_params.initial_lapse_type <= 4),
                        "must be 0, 1, 2, 3, or 4");

        check_parameter("wormhole_throat_radius", wormhole_params.b0,
                        wormhole_params.b0 > 0.0,
                        "must be positive");

        check_parameter("wormhole_phi_perturbation_width",
                        wormhole_params.phi_perturbation_width,
                        (wormhole_params.phi_perturbation_amplitude == 0.0 &&
                         wormhole_params.phi_monopole_amplitude == 0.0) ||
                            (wormhole_params.phi_perturbation_width > 0.0),
                        "must be > 0 when any phi perturbation amplitude is nonzero");

        check_parameter("wormhole_matter_model", wormhole_matter_model,
                        wormhole_matter_model == "exotic_scalar" ||
                            wormhole_matter_model == "no_matter" ||
                            wormhole_matter_model == "effective_teo" ||
                            wormhole_matter_model == "dust" ||
                            wormhole_matter_model == "oscillon_scalar" ||
                            wormhole_matter_model == "complex_scalar",
                        "must be exotic_scalar, no_matter, effective_teo, dust, "
                        "oscillon_scalar, or complex_scalar");
    }

    bool calculate_constraint_norms{};
    std::string wormhole_matter_model;

    std::string recipe_initial_data_file;
    ExternalGridInitialData::params_t external_grid_params{};

    SupportedWormholeInitialData::params_t wormhole_params{};

    // Trajectory-guided matter pump (Rung 2 active support).
    int trajectory_mode{0};
    TrajectoryParams trajectory_params{};
    double trajectory_pump_frequency{0.0};
    // Per-lump action state populated from trajectory params each coarse step.
    std::array<double, RL_MAX_LUMPS> rl_pump_amplitude{};
    std::array<double, RL_MAX_LUMPS> rl_pump_frequency{};
    std::array<double, RL_MAX_LUMPS> rl_pump_phase{};
    double rl_pump_width{1.5};
    double rl_pump_max_amplitude{0.05};
    double rl_l2_ham_governor_center{0.035};
    double rl_l2_ham_governor_width{0.003};
    double rl_pump_kp{0.0};
    double rl_pump_kd{0.0};
    int rl_pump_target_profile{0};
    double rl_pump_target_width{0.0};
    double rl_pump_target_amp{0.0};
    // Pump ramp schedule (Phase 6 collapse trigger).  t_start < 0 => never.
    double pump_ramp_t_start{-1.0};
    double pump_ramp_t_end{0.0};
    double pump_ramp_floor{0.0};

    // Support-strength ramp schedule (Phase 6 collapse-on-command trigger).
    double support_ramp_t_start{-1.0};
    double support_ramp_t_end{0.0};
    double support_ramp_floor{0.0};

    // Numerical sponge zone (radially-ramped extra KO dissipation).
    SpongeZoneParams sponge_params{};

  private:
    void load_trajectory_lump(GRParmParse &pp, int k, PerLumpTrajectory &lk)
    {
        std::ostringstream pfx;
        pfx << "trajectory_lump" << k << "_";
        const std::string p = pfx.str();
        pp.load((p + "R0").c_str(), lk.R0, 5.0);
        pp.load((p + "omega_rot").c_str(), lk.omega_rot, 0.0);
        pp.load((p + "phase0").c_str(), lk.phase0, 0.0);
        pp.load((p + "tilt_theta").c_str(), lk.tilt_theta, 0.0);
        pp.load((p + "tilt_phi").c_str(), lk.tilt_phi, 0.0);
        pp.load((p + "well_depth").c_str(), lk.well_depth, 0.05);
        pp.load((p + "v_rad").c_str(), lk.v_rad, 0.0);
    }
};

#endif /* SIMULATIONPARAMETERS_HPP */