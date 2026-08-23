#ifndef SIMULATIONPARAMETERS_HPP
#define SIMULATIONPARAMETERS_HPP

#include "ExternalGridInitialData.hpp"
#include "GRParmParse.hpp"
#include "GRTresnaScalarLayout.hpp"
#include "GridTreadmill.hpp"
#include "RadialRecipeInitialData.hpp"
#include "SimulationParametersBase.hpp"
#include "SpongeZone.hpp"
#include "TrajectoryParams.hpp"

#include <AMReX.H>

#include <array>
#include <sstream>
#include <string>

class SimulationParameters : public SimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        read_shared_params(pp);
        read_recipe_params(pp);
        check_params();
    }

    void read_shared_params(GRParmParse &pp)
    {
        pp.load("calculate_constraint_norms", calculate_constraint_norms, false);
        pp.load("calculate_energy_conditions", calculate_energy_conditions,
                false);
        pp.load("calculate_curvature_invariants",
                calculate_curvature_invariants, false);

        // Matter sector selection.  The constrained initial-data recipe may
        // require exotic (phantom, rho <= 0) matter to source a wormhole/warp
        // geometry.  When recipe_exotic_matter is set the level evolves an
        // ExoticScalarField whose stress-energy is the canonical one scaled by
        // -recipe_support_strength, so the matter that is actually evolved
        // matches the matter the geometry was reconstructed for.
        pp.load("recipe_exotic_matter", recipe_exotic_matter, false);
        pp.load("recipe_support_strength", recipe_support_strength, 1.0);

        pp.load("recipe_matter_model", recipe_matter_model, std::string(""));
        pp.load("recipe_num_scalar_fields", recipe_num_scalar_fields, 0);
        pp.load("recipe_scalar_mass", recipe_scalar_mass, 0.0);
        pp.load("recipe_scalar_lambda", recipe_scalar_lambda, 0.0);
        pp.load("recipe_scalar_mu", recipe_scalar_mu, 0.0);
        pp.load("recipe_scalar_sign", recipe_scalar_sign, 1.0);
        load_scalar_field_signs(pp);

        pp.load("recipe_initial_data_file", recipe_initial_data_file,
                std::string(""));
        if (!recipe_initial_data_file.empty())
        {
            external_grid_params.gridinit_file = recipe_initial_data_file;
            external_grid_params.grid_center = center;
        }

        pp.load("rl_enabled", rl_enabled, false);
        pp.load("rl_coarse_step_interval", rl_coarse_step_interval, 64);
        pp.load("rl_zmq_port", rl_zmq_port, 5555);
        pp.load("rl_zmq_timeout_ms", rl_zmq_timeout_ms, 120000);
        pp.load("rl_num_lumps", rl_num_lumps, 1);
        if (rl_num_lumps < 1)
            rl_num_lumps = 1;
        // Refuse rather than truncate -- same reasoning as trajectory_num_lumps
        // below (DebugPreGPU.md PG-9): a silently reduced lump count makes every
        // rl_lump<k>_* key past the cap a dead knob that still looks configured.
        if (rl_num_lumps > GRTRESNA_MAX_INDEPENDENT_SCALARS)
        {
            amrex::Abort(
                "rl_num_lumps = " + std::to_string(rl_num_lumps) +
                " exceeds GRTRESNA_MAX_INDEPENDENT_SCALARS = " +
                std::to_string(GRTRESNA_MAX_INDEPENDENT_SCALARS) +
                ". Lower rl_num_lumps or raise the cap in "
                "Source/Matter/GRTresnaScalarLayout.hpp.");
        }
        pp.load("rl_pump_width", rl_pump_width, 1.5);
        pp.load("rl_pump_max_amplitude", rl_pump_max_amplitude, 0.05);
        pp.load("rl_l2_ham_governor_center", rl_l2_ham_governor_center, 0.035);
        pp.load("rl_l2_ham_governor_width", rl_l2_ham_governor_width, 0.003);
        // Closed-loop PD "trap" controller gains for the matter pump.  When
        // k_p > 0 the pump drives the field toward the moving target soliton
        // (error-proportional, self-limiting) instead of the legacy open-loop
        // source.  Default 0 => legacy source pump (backward compatible).
        pp.load("rl_pump_kp", rl_pump_kp, 0.0);
        pp.load("rl_pump_kd", rl_pump_kd, 0.0);
        // Per-sector gain override for the phantom (Phi-) field.  < 0 => the
        // phantom sector inherits the shared gains above (bit-identical to the
        // historical single-gain behaviour).  See RLMatterPumpParams for the
        // measured decay rates that motivate a separate phantom gain.
        pp.load("rl_pump_kp_phantom", rl_pump_kp_phantom, -1.0);
        pp.load("rl_pump_kd_phantom", rl_pump_kd_phantom, -1.0);
        // Trap-target matter shape: 0 Gaussian (legacy), 2 sech bound lump.
        // target width = physical bound size 1/sqrt(m^2-omega^2); <=0 => pump width.
        pp.load("rl_pump_target_profile", rl_pump_target_profile, 0);
        pp.load("rl_pump_target_width", rl_pump_target_width, 0.0);
        // Trap-target central amplitude (self-grav boson star: bs_phi_c).  0 =>
        // use the per-site amplitude (legacy).
        pp.load("rl_pump_target_amp", rl_pump_target_amp, 0.0);
        // Superposed-target PD law: overlapping same-sector sites share one
        // summed target and a capped weight instead of summing per-site
        // errors (which strips matter where lumps overlap).  0 = legacy.
        pp.load("rl_pump_superpose_targets", rl_pump_superpose_targets, 0);
        // Transient "igniter" pump: when >= 0 and time >= stop, force the pump
        // fully off (k_p=k_d=0, num_sites=0) so post-stop evolution is
        // conservative Einstein-Klein-Gordon.  Default -1 => never auto-stop.
        pp.load("rl_pump_stop_time", rl_pump_stop_time, -1.0);
        // Controller energy-momentum reservoir:
        //   0 = off (default; bit-identical to pre-reservoir runs)
        //   1 = ledger (evolve reservoir; include in constraint diagnostic)
        //   2 = backreaction (also include reservoir in the CCZ4 RHS)
        pp.load("controller_reservoir_mode", controller_reservoir_mode, 0);
        if (controller_reservoir_mode < 0)
            controller_reservoir_mode = 0;
        if (controller_reservoir_mode > 2)
            controller_reservoir_mode = 2;
        // Per-lump action state starts at zero; the RL agent populates it.
        rl_pump_amplitude.fill(0.0);
        rl_pump_frequency.fill(0.0);
        rl_pump_phase.fill(0.0);
        // Seed lump centres from flat "x0 y0 z0 x1 y1 z1 ..." lists (elite geom).
        load_rl_lump_seed_axis(pp, "rl_lump_seed_x", rl_lump_seed_x);
        load_rl_lump_seed_axis(pp, "rl_lump_seed_y", rl_lump_seed_y);
        load_rl_lump_seed_axis(pp, "rl_lump_seed_z", rl_lump_seed_z);

        // Trajectory-guided geometry survey (opt-in; trajectory_mode=0 is off).
        // When active, the pump site centres are overridden each coarse step by
        // parametric trajectory equations evaluated on the CPU (no ZMQ, no RL).
        pp.load("trajectory_mode", trajectory_mode, 0);
        {
            int n_traj = 5;
            pp.load("trajectory_num_lumps", n_traj, 5);
            if (n_traj < 1)
                n_traj = 1;
            // Refuse rather than truncate.  This used to silently clamp to the
            // compile-time cap, so a campaign configured for 8 lumps drove 5 and
            // said nothing: trajectory_lump5..7_* were written into every
            // params.txt and read by no one, leaving 21 declared search
            // dimensions with zero effect on the evolution while the archive
            // reported coverage over all of them.  A run that cannot do what its
            // config asks must not start.  See DebugPreGPU.md PG-9.
            if (n_traj > GRTRESNA_MAX_INDEPENDENT_SCALARS)
            {
                amrex::Abort(
                    "trajectory_num_lumps = " + std::to_string(n_traj) +
                    " exceeds GRTRESNA_MAX_INDEPENDENT_SCALARS = " +
                    std::to_string(GRTRESNA_MAX_INDEPENDENT_SCALARS) +
                    ". The extra lumps' trajectory parameters would be ignored "
                    "silently. Either lower trajectory_num_lumps or raise the "
                    "cap in Source/Matter/GRTresnaScalarLayout.hpp (which also "
                    "widens the state vector by 2 components per lump).");
            }
            trajectory_params.num_lumps = n_traj;
        }

        // Shared trajectory parameters.
        pp.load("trajectory_A_breath", trajectory_params.A_breath, 0.0);
        pp.load("trajectory_omega_breath", trajectory_params.omega_breath,
                0.0);
        pp.load("trajectory_z_amp", trajectory_params.z_amp, 0.0);
        pp.load("trajectory_omega_z", trajectory_params.omega_z, 0.0);
        pp.load("trajectory_well_width", trajectory_params.well_width, 1.5);
        // Shared pump frequency for complex-scalar trajectory mode.
        // When nonzero the pump drives Pi and Pi2 with this angular frequency
        // (should match bs_omega for coherent U(1) charge injection).
        pp.load("trajectory_pump_frequency", trajectory_pump_frequency, 0.0);

        // Per-lump trajectory parameters: trajectory_lump{k}_*.
        for (int k = 0; k < trajectory_params.num_lumps; ++k)
        {
            auto &lk = trajectory_params.lumps[k];
            load_trajectory_lump(pp, k, lk);
        }

        // Numerical sponge zone: radially-ramped extra KO dissipation in an
        // outer shell to suppress boundary reflections without enlarging the box.
        pp.load("sponge_enabled", sponge_params.enabled, false);
        pp.load("sponge_inner_radius", sponge_params.inner_radius, 24.0);
        pp.load("sponge_outer_radius", sponge_params.outer_radius, 32.0);
        pp.load("sponge_strength", sponge_params.strength, 4.0);
        pp.load("sponge_ramp_power", sponge_params.ramp_power, 4);
        pp.load("sponge_center", sponge_params.center, center);

        // Recentring box ("treadmill"): once the source has drifted, carry the
        // DATA back toward the box centre by a whole number of cells and keep
        // an odometer of how far, so a runaway can be followed for as long as
        // wanted without enlarging the box.  Off by default, so every archived
        // cell is bit-for-bit unaffected.  Design and validation ladder:
        // research/bondi_dipole/docs/CHASE_TO_03C.md.
        pp.load("treadmill_enabled", treadmill_params.enabled, false);
        pp.load("treadmill_axis", treadmill_params.axis, 0);
        pp.load("treadmill_threshold", treadmill_params.threshold, 2.0);
        pp.load("treadmill_check_interval", treadmill_params.check_interval,
                20);
        pp.load("treadmill_ball_radius", treadmill_params.ball_radius, 8.0);
        pp.load("treadmill_fill_mode", treadmill_params.fill_mode, 0);
    }

    void read_recipe_params(GRParmParse &pp)
    {
        pp.load("recipe_num_bases", recipe_params.num_bases, 4);
        pp.load("center", recipe_params.grid_center, center);
        pp.load("recipe_basis_width", recipe_params.basis_width, 1.0);
        pp.load("recipe_basis_radius_max", recipe_params.basis_radius_max, 8.0);

        pp.load("recipe_chi_asymptotic", recipe_params.chi_asymptotic, 1.0);
        pp.load("recipe_alpha_asymptotic", recipe_params.alpha_asymptotic, 1.0);
        pp.load("recipe_beta_asymptotic", recipe_params.beta_asymptotic, 0.0);
        pp.load("recipe_K_asymptotic", recipe_params.K_asymptotic, 0.0);
        pp.load("recipe_phi_asymptotic", recipe_params.phi_asymptotic, 0.0);
        pp.load("recipe_Pi_asymptotic", recipe_params.Pi_asymptotic, 0.0);

        load_coeff_array(pp, "recipe_chi_coeff", recipe_params.chi_coeffs);
        load_coeff_array(pp, "recipe_alpha_coeff", recipe_params.alpha_coeffs);
        load_coeff_array(pp, "recipe_beta_coeff", recipe_params.beta_coeffs);
        load_coeff_array(pp, "recipe_K_coeff", recipe_params.K_coeffs);
        load_coeff_array(pp, "recipe_phi_coeff", recipe_params.phi_coeffs);
        load_coeff_array(pp, "recipe_Pi_coeff", recipe_params.Pi_coeffs);

        pp.load("recipe_num_chi_angular_modes",
                recipe_params.num_chi_angular_modes, 0);
        load_angular_modes(pp, "recipe_chi_mode",
                           recipe_params.num_chi_angular_modes,
                           recipe_params.chi_angular_modes);

        pp.load("recipe_num_chi_Ylm_modes", recipe_params.num_chi_Ylm_modes, 0);
        load_Ylm_modes(pp);

        pp.load("recipe_num_lapse_angular_modes",
                recipe_params.num_lapse_angular_modes, 0);
        load_angular_modes(pp, "recipe_lapse_mode",
                           recipe_params.num_lapse_angular_modes,
                           recipe_params.lapse_angular_modes);

        pp.load("recipe_num_beta_angular_modes",
                recipe_params.num_beta_angular_modes, 0);
        load_angular_modes(pp, "recipe_beta_mode",
                           recipe_params.num_beta_angular_modes,
                           recipe_params.beta_angular_modes);

        pp.load("recipe_num_K_angular_modes",
                recipe_params.num_K_angular_modes, 0);
        load_angular_modes(pp, "recipe_K_mode",
                           recipe_params.num_K_angular_modes,
                           recipe_params.K_angular_modes);
    }

    void check_params()
    {
        check_parameter("recipe_num_bases", recipe_params.num_bases,
                        recipe_params.num_bases > 0 &&
                            recipe_params.num_bases <=
                                RadialRecipeInitialData::MAX_BASES,
                        "must be between 1 and MAX_BASES");

        check_parameter("recipe_basis_width", recipe_params.basis_width,
                        recipe_params.basis_width > 0.0, "must be positive");

        check_parameter("recipe_basis_radius_max",
                        recipe_params.basis_radius_max,
                        recipe_params.basis_radius_max > 0.0,
                        "must be positive");

        check_parameter("recipe_num_chi_angular_modes",
                        recipe_params.num_chi_angular_modes,
                        recipe_params.num_chi_angular_modes >= 0 &&
                            recipe_params.num_chi_angular_modes <=
                                RadialRecipeInitialData::MAX_ANGULAR_MODES,
                        "must be between 0 and MAX_ANGULAR_MODES");

        check_parameter("recipe_num_chi_Ylm_modes",
                        recipe_params.num_chi_Ylm_modes,
                        recipe_params.num_chi_Ylm_modes >= 0 &&
                            recipe_params.num_chi_Ylm_modes <=
                                RadialRecipeInitialData::MAX_YLM_MODES,
                        "must be between 0 and MAX_YLM_MODES");

        check_angular_mode_count("recipe_num_lapse_angular_modes",
                                 recipe_params.num_lapse_angular_modes);
        check_angular_mode_count("recipe_num_beta_angular_modes",
                                 recipe_params.num_beta_angular_modes);
        check_angular_mode_count("recipe_num_K_angular_modes",
                                 recipe_params.num_K_angular_modes);
    }

    void check_angular_mode_count(const char *name, int count)
    {
        check_parameter(name, count,
                        count >= 0 &&
                            count <= RadialRecipeInitialData::MAX_ANGULAR_MODES,
                        "must be between 0 and MAX_ANGULAR_MODES");
    }

    bool calculate_constraint_norms{};
    bool calculate_energy_conditions{};
    bool calculate_curvature_invariants{};
    bool recipe_exotic_matter{};
    double recipe_support_strength{1.0};

    std::string recipe_matter_model;
    int recipe_num_scalar_fields{};
    std::array<int, GRTRESNA_MAX_INDEPENDENT_SCALARS> recipe_scalar_field_signs{};
    double recipe_scalar_mass{};
    double recipe_scalar_lambda{};
    double recipe_scalar_mu{};
    double recipe_scalar_sign{1.0};

    std::string recipe_initial_data_file;
    ExternalGridInitialData::params_t external_grid_params{};

    RadialRecipeInitialData::params_t recipe_params{};

    // RL closed-loop control (opt-in; defaults preserve IVP behaviour)
    bool rl_enabled{};
    int rl_coarse_step_interval{64};
    int rl_zmq_port{5555};
    int rl_zmq_timeout_ms{30000};
    int rl_num_lumps{1};
    // Per-lump action state (set by the RL agent at runtime, persisted across
    // intervals).  Index k corresponds to lump k.
    std::array<double, GRTRESNA_MAX_INDEPENDENT_SCALARS> rl_pump_amplitude{};
    std::array<double, GRTRESNA_MAX_INDEPENDENT_SCALARS> rl_pump_frequency{};
    std::array<double, GRTRESNA_MAX_INDEPENDENT_SCALARS> rl_pump_phase{};
    // Initial lump centres (tracker search-ball seeds) from elite geometry.
    std::array<double, GRTRESNA_MAX_INDEPENDENT_SCALARS> rl_lump_seed_x{};
    std::array<double, GRTRESNA_MAX_INDEPENDENT_SCALARS> rl_lump_seed_y{};
    std::array<double, GRTRESNA_MAX_INDEPENDENT_SCALARS> rl_lump_seed_z{};
    double rl_pump_width{1.5};
    double rl_pump_max_amplitude{0.05};
    double rl_l2_ham_governor_center{0.035};
    double rl_l2_ham_governor_width{0.003};
    double rl_pump_kp{0.0}; //!< PD trap proportional gain (0 => legacy source)
    double rl_pump_kd{0.0}; //!< PD trap derivative gain
    double rl_pump_kp_phantom{-1.0}; //!< phantom-sector k_p (<0 => inherit k_p)
    double rl_pump_kd_phantom{-1.0}; //!< phantom-sector k_d (<0 => inherit k_d)
    int rl_pump_target_profile{0}; //!< 0 Gaussian, 2 sech bound-lump target
    double rl_pump_target_width{0.0}; //!< physical 1/sqrt(m^2-omega^2); <=0 => width
    double rl_pump_target_amp{0.0}; //!< trap-target central amplitude (bs_phi_c); <=0 => site amp
    int rl_pump_superpose_targets{0}; //!< 1 => summed sector target + capped weight; 0 legacy
    double rl_pump_stop_time{-1.0}; //!< >=0 => pump off for time>=stop; -1 => never
    int controller_reservoir_mode{0}; //!< 0 off, 1 ledger, 2 backreaction

    // Trajectory-guided geometry survey (Independent of RL; no ZMQ needed).
    int trajectory_mode{0};        //!< 0 = off, 1 = parametric trajectory
    double trajectory_pump_frequency{0.0}; //!< Shared pump angular frequency (bs_omega for complex scalar)
    TrajectoryParams trajectory_params{};

    // Numerical sponge zone (radially-ramped extra KO dissipation).
    SpongeZoneParams sponge_params{};

    // Recentring box (exact whole-cell translation + odometer).
    GridTreadmillParams treadmill_params{};

  private:
    // NOTE on multi-value keys: `key = 1 -1 -1 1 -1` is tokenized by
    // ParmParse into five values, and a scalar query returns ONLY token 0.
    // Both helpers below used to read the key into ONE std::string and split
    // it -- which silently kept just the first value. For
    // recipe_scalar_field_signs that parsed `1 -1 -1 1 -1` as [1, 0, 0, 0, 0],
    // and since 0 routes to canonical, EVERY pump spotlight drove the
    // canonical field: the phantom sector never received any pump force on
    // any bicomplex campaign up to 2026-07-28 (Debug.md 19.8). Multi-value
    // keys must go through countval + getarr.
    void load_rl_lump_seed_axis(
        GRParmParse &pp, const char *key,
        std::array<double, GRTRESNA_MAX_INDEPENDENT_SCALARS> &out)
    {
        out.fill(0.0);
        const int nvals = pp.countval(key);
        if (nvals < 1)
            return;
        std::vector<double> vals;
        pp.getarr(key, vals, 0,
                  std::min(nvals, GRTRESNA_MAX_INDEPENDENT_SCALARS));
        for (std::size_t i = 0; i < vals.size(); ++i)
        {
            out[i] = vals[i];
        }
    }

    void load_scalar_field_signs(GRParmParse &pp)
    {
        const int nvals = pp.countval("recipe_scalar_field_signs");
        if (nvals > 0)
        {
            std::vector<int> signs;
            pp.getarr("recipe_scalar_field_signs", signs, 0,
                      std::min(nvals, GRTRESNA_MAX_INDEPENDENT_SCALARS));
            for (std::size_t i = 0; i < signs.size(); ++i)
            {
                recipe_scalar_field_signs[i] = signs[i];
            }
        }
        else
        {
            for (int k = 0; k < GRTRESNA_MAX_INDEPENDENT_SCALARS; ++k)
            {
                std::ostringstream key;
                key << "recipe_scalar_field_sign_" << k;
                pp.load(key.str().c_str(), recipe_scalar_field_signs[k], 1);
            }
        }
        // Echo what was actually parsed: the silent first-token bug above
        // survived four campaigns because nothing ever printed this.
        amrex::Print() << "recipe_scalar_field_signs parsed:";
        for (int k = 0; k < GRTRESNA_MAX_INDEPENDENT_SCALARS; ++k)
        {
            amrex::Print() << " " << recipe_scalar_field_signs[k];
        }
        amrex::Print() << "\n";
    }

    void load_coeff_array(GRParmParse &pp, const char *prefix,
                          std::array<double, RadialRecipeInitialData::MAX_BASES>
                              &coeffs)
    {
        for (int n = 0; n < recipe_params.num_bases; ++n)
        {
            std::ostringstream key;
            key << prefix << "_" << n;
            pp.load(key.str().c_str(), coeffs[n], 0.0);
        }
    }

    void load_angular_modes(
        GRParmParse &pp, const char *prefix, int num_modes,
        std::array<RadialRecipeInitialData::AngularMode,
                   RadialRecipeInitialData::MAX_ANGULAR_MODES> &modes)
    {
        for (int n = 0; n < num_modes; ++n)
        {
            auto &mode = modes[n];
            std::ostringstream ell_key, amp_key, rc_key, rw_key;
            ell_key << prefix << "_ell_" << n;
            amp_key << prefix << "_amp_" << n;
            rc_key << prefix << "_rc_" << n;
            rw_key << prefix << "_rw_" << n;
            pp.load(ell_key.str().c_str(), mode.ell, 0);
            pp.load(amp_key.str().c_str(), mode.amplitude, 0.0);
            pp.load(rc_key.str().c_str(), mode.radial_center, 0.0);
            pp.load(rw_key.str().c_str(), mode.radial_width, 1.0);
        }
    }

    void load_Ylm_modes(GRParmParse &pp)
    {
        for (int n = 0; n < recipe_params.num_chi_Ylm_modes; ++n)
        {
            auto &mode = recipe_params.chi_Ylm_modes[n];
            std::ostringstream l_key, m_key, amp_key, rc_key, rw_key;
            l_key << "recipe_chi_Ylm_l_" << n;
            m_key << "recipe_chi_Ylm_m_" << n;
            amp_key << "recipe_chi_Ylm_amp_" << n;
            rc_key << "recipe_chi_Ylm_rc_" << n;
            rw_key << "recipe_chi_Ylm_rw_" << n;
            pp.load(l_key.str().c_str(), mode.ell, 0);
            pp.load(m_key.str().c_str(), mode.em, 0);
            pp.load(amp_key.str().c_str(), mode.amplitude, 0.0);
            pp.load(rc_key.str().c_str(), mode.radial_center, 0.0);
            pp.load(rw_key.str().c_str(), mode.radial_width, 1.0);
        }
    }

    //! Load per-lump trajectory: trajectory_lump{k}_R0, _omega_rot, etc.
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
