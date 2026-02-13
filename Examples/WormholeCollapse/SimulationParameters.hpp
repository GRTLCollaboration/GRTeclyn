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
        // Read Wormhole specific parameters
        pp.load("wormhole_throat_radius", wormhole_params.throat_radius, 1.0);
        pp.load("wormhole_k_amplitude", wormhole_params.k_amplitude, 0.1);
        pp.load("wormhole_k_width", wormhole_params.k_width, 1.0);

        // Center of the wormhole
        pp.load("center", wormhole_params.center, center);
    }

    void check_params()
    {
        check_parameter("wormhole_throat_radius", wormhole_params.throat_radius,
                        wormhole_params.throat_radius > 0.0,
                        "must be positive");
    }

    bool calculate_constraint_norms{};

    // Store parameters for Initial Data
    WormholeInitialData::params_t wormhole_params{};
};

#endif /* SIMULATIONPARAMETERS_HPP */