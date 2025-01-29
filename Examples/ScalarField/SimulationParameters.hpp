/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"

// Problem specific includes:
#include "InitialBackgroundData.hpp"
#include "InitialScalarData.hpp"
#include "Potential.hpp"
#include "RandomField.hpp"

class SimulationParameters : public SimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        // read the problem specific params
        read_params(pp);
        check_params();
    }

    void read_params(GRParmParse &pp)
    {
	    initial_params.center =
            center; // already read in SimulationParametersBase
         pp.load("G_Newton", G_Newton,
                 0.0); // for now the example neglects backreaction
	    pp.load("scalar_mass", potential_params.scalar_mass, 0.1);	

	    pp.load("scalar_amplitude", background_params.phi0, 0.0);
	    pp.load("scalar_velocity", background_params.Pi0, 0.0);
        pp.load("scalar_mass", background_params.m, 0.0);

        pp.load("num_scalar_fields", random_field_params.num_scalar_fields, 0);
        pp.load("calc_tensor_field", random_field_params.calc_tensor_field, 0);
        pp.load("L_full", random_field_params.L, 1.);
        pp.load("A", random_field_params.A, 1.);
        pp.load("N_full", random_field_params.N_readin, 32);
        pp.load("N_fine", random_field_params.N_fine, random_field_params.N_readin);
        pp.load("use_rand", random_field_params.use_rand, 1);
        pp.load("random_seed", random_field_params.random_seed, 3539263);

        pp.load("use_window", random_field_params.use_window, 0);
        pp.load("kstar", random_field_params.kstar, 0.);
        pp.load("Delta", random_field_params.Delta, 1.);

        pp.load("calc_binned_power_spectrum", random_field_params.calc_binned_power_spectrum, 0);
        pp.load("bin_number", random_field_params.bin_number, random_field_params.N_readin/2); 
        pp.load("calc_config_space_mode_fns", random_field_params.calc_config_space_mode_fns, 0);
        pp.load("calc_higher_order_statistics", random_field_params.calc_higher_order_statistics, 0);
    }

    void check_params()
    {
        warn_parameter("scalar_mass", background_params.m,
                       background_params.m <
                           0.2 / coarsest_dx / dt_multiplier,
                       "oscillations of scalar field do not appear to be "
                       "resolved on coarsest level");

        warn_parameter("kstar", random_field_params.kstar,
                       random_field_params.kstar > 0,
                       "cut-off frequency index must be positive");

        check_parameter("Delta", random_field_params.Delta,
                       random_field_params.Delta > 0,
                       "cut-off width must be positive and non-zero");
    }

    // Initial data for matter and potential and BH
    double G_Newton;
    Potential::params_t potential_params;
    InitialBackgroundData::params_t background_params;
    InitialScalarData::params_t initial_params;
    RandomField::params_t random_field_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
