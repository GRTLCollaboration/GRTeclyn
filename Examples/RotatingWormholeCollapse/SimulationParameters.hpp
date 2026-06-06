#ifndef SIMULATIONPARAMETERS_HPP
#define SIMULATIONPARAMETERS_HPP

#include "ExternalGridInitialData.hpp"
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

        pp.load("wormhole_support_strength", wormhole_params.support_strength,
                1.0);

        pp.load("wormhole_phi_monopole_amplitude",
                wormhole_params.phi_monopole_amplitude, 0.0);
        pp.load("wormhole_phi_perturbation_amplitude",
                wormhole_params.phi_perturbation_amplitude, 0.0);
        pp.load("wormhole_phi_perturbation_width",
                wormhole_params.phi_perturbation_width, 0.0);

        pp.load("recipe_initial_data_file", recipe_initial_data_file,
                std::string(""));
        if (!recipe_initial_data_file.empty())
        {
            external_grid_params.gridinit_file = recipe_initial_data_file;
            external_grid_params.grid_center = center;
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
                            wormhole_matter_model == "oscillon_scalar",
                        "must be exotic_scalar, no_matter, effective_teo, dust, or oscillon_scalar");
    }

    bool calculate_constraint_norms{};
    std::string wormhole_matter_model;

    std::string recipe_initial_data_file;
    ExternalGridInitialData::params_t external_grid_params{};

    SupportedWormholeInitialData::params_t wormhole_params{};
};

#endif /* SIMULATIONPARAMETERS_HPP */