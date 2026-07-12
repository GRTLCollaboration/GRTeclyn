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
};

#endif /* SIMULATIONPARAMETERS_HPP */