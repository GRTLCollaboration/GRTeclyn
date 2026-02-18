#ifndef SIMULATIONPARAMETERS_HPP
#define SIMULATIONPARAMETERS_HPP

// General includes
#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"
#include "WormholeInitialData.hpp" // Add this

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
        // Grid center for coordinate mapping
        pp.load("center", wormhole_params.grid_center, center);

        // Backward-compatible single value
        double b0_single = 1.0;
        std::array<double, AMREX_SPACEDIM> c_single = {0.0, 0.0, 0.0};
        pp.load("wormhole_throat_radius", b0_single, 1.0);
        pp.load("wormhole_center", c_single,
                std::array<double, AMREX_SPACEDIM>{2.0, 0.0, 0.0});

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
    }

    void check_params()
    {
        check_parameter("wormhole_throat_radius_A", wormhole_params.throat_radius_A,
                        wormhole_params.throat_radius_A > 0.0,
                        "must be positive");
        check_parameter("wormhole_throat_radius_B", wormhole_params.throat_radius_B,
                        wormhole_params.throat_radius_B > 0.0,
                        "must be positive");
    }

    bool calculate_constraint_norms{};

    // Store parameters for Initial Data
    WormholeInitialData::params_t wormhole_params{};
};

#endif /* SIMULATIONPARAMETERS_HPP */