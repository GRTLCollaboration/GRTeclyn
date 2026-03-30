#ifndef SIMULATIONPARAMETERS_HPP
#define SIMULATIONPARAMETERS_HPP

#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"
#include "SupportedWormholeInitialData.hpp"

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

        pp.load("wormhole_support_strength", wormhole_params.support_strength,
                1.0);

        pp.load("wormhole_phi_monopole_amplitude",
                wormhole_params.phi_monopole_amplitude, 0.0);
        pp.load("wormhole_phi_perturbation_amplitude",
                wormhole_params.phi_perturbation_amplitude, 0.0);
        pp.load("wormhole_phi_perturbation_width",
                wormhole_params.phi_perturbation_width, 0.0);
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
    }

    bool calculate_constraint_norms{};

    SupportedWormholeInitialData::params_t wormhole_params{};
};

#endif /* SIMULATIONPARAMETERS_HPP */