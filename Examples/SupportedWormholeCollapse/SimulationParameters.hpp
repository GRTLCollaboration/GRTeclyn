#ifndef SIMULATIONPARAMETERS_HPP
#define SIMULATIONPARAMETERS_HPP

// General includes
#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"
#include "SupportedWormholeInitialData.hpp" // Add this

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
        // Keep constraint calculation if needed
        pp.load("calculate_constraint_norms", calculate_constraint_norms,
                false);
    }

    void read_wormhole_params(GRParmParse &pp)
    {
        // Select between different initial-data realisations
        pp.load("wormhole_metric_type", wormhole_params.metric_type, 0);
        pp.load("wormhole_initial_lapse_type", wormhole_params.initial_lapse_type,
                0);

        // Grid center for coordinate mapping
        pp.load("center", wormhole_params.grid_center, center);

        // Backward-compatible single value
        double b0_single = 1.0;
        // Default to a two-mouth separation of 30 (centred about the origin in
        // physical coordinates when `center = L_full/2`).
        std::array<double, AMREX_SPACEDIM> c_single = {15.0, 0.0, 0.0};
        pp.load("wormhole_throat_radius", b0_single, 1.0);
        pp.load("wormhole_center", c_single,
                std::array<double, AMREX_SPACEDIM>{15.0, 0.0, 0.0});

        // Two-mouth parameters (preferred)
        pp.load("wormhole_throat_radius_A", wormhole_params.throat_radius_A,
                b0_single);
        pp.load("wormhole_throat_radius_B", wormhole_params.throat_radius_B,
                b0_single);
        pp.load("wormhole_centerA", wormhole_params.centerA, c_single);

        std::array<double, AMREX_SPACEDIM> default_centerB = {
            -wormhole_params.centerA[0],
            -wormhole_params.centerA[1],
            -wormhole_params.centerA[2],
        };
        pp.load("wormhole_centerB", wormhole_params.centerB, default_centerB);

        // Legacy/debug option
        pp.load("wormhole_use_cartesian_gamma", wormhole_params.use_cartesian_gamma,
                false);

        pp.load("wormhole_support_strength", wormhole_params.support_strength,
                1.0);
        pp.load("wormhole_support_ramp_start", wormhole_params.support_ramp_start,
                -1.0); // negative means no ramp
        pp.load("wormhole_support_ramp_duration", wormhole_params.support_ramp_duration,
                5.0);

        // Kick parameters:
        // - wormhole_k_monopole_amplitude controls the coherent compressive part
        // - wormhole_k_quadrupole_amplitude controls the Y20 deformation
        // - wormhole_k_amplitude is kept as a backward-compatible alias for the
        //   quadrupole piece
        pp.load("wormhole_k_amplitude", wormhole_params.k_amplitude, 0.0);
        pp.load("wormhole_k_monopole_amplitude",
                wormhole_params.k_monopole_amplitude, 0.0);
        pp.load("wormhole_k_quadrupole_amplitude",
                wormhole_params.k_quadrupole_amplitude,
                wormhole_params.k_amplitude);
        pp.load("wormhole_k_width", wormhole_params.k_width, 0.0);
    }

    void check_params()
    {
        check_parameter("wormhole_metric_type", wormhole_params.metric_type,
                        (wormhole_params.metric_type == 0) ||
                            (wormhole_params.metric_type == 1),
                        "must be 0 (two-mouth) or 1 (single-throat)");
        check_parameter("wormhole_initial_lapse_type",
                        wormhole_params.initial_lapse_type,
                        (wormhole_params.initial_lapse_type >= 0) &&
                            (wormhole_params.initial_lapse_type <= 2),
                        "must be 0, 1, or 2");

        check_parameter("wormhole_throat_radius_A", wormhole_params.throat_radius_A,
                        wormhole_params.throat_radius_A > 0.0,
                        "must be positive");
        check_parameter("wormhole_throat_radius_B", wormhole_params.throat_radius_B,
                        wormhole_params.throat_radius_B > 0.0,
                        "must be positive");

        // Kick params: width must be > 0 if amplitude is nonzero
        check_parameter("wormhole_k_width", wormhole_params.k_width,
                        ((wormhole_params.k_monopole_amplitude == 0.0) &&
                         (wormhole_params.k_quadrupole_amplitude == 0.0)) ||
                            (wormhole_params.k_width > 0.0),
                        "must be > 0 when any wormhole kick amplitude is nonzero");
    }

    bool calculate_constraint_norms{};

    // Store parameters for Initial Data
    SupportedWormholeInitialData::params_t wormhole_params{};
};

#endif /* SIMULATIONPARAMETERS_HPP */