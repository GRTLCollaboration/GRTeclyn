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
	    pp.load("potential_param_1", potential_params.scalar_mass, 0.1);	

	    pp.load("background_phi", background_params.phi0, 0.0);
	    pp.load("background_dphi", background_params.Pi0, 0.0);
        pp.load("potential_param_1", background_params.m, 0.0);
        pp.load("G_Newton", background_params.G_Newton, 1.);

        pp.load("tensor_init", random_field_params.tensor_init, 1);
        pp.load("L", random_field_params.L, 1.);
        pp.load("A", random_field_params.A, 1.);
        pp.load("N", random_field_params.N_readin, 32);
        //pp.load("N_fine", random_field_params.N_fine, random_field_params.N_readin);
        pp.load("use_rand", random_field_params.use_rand, 1);
        pp.load("random_seed", random_field_params.random_seed, 3539263);

        pp.load("use_window", random_field_params.use_window, 0);
        pp.load("kstar", random_field_params.kstar, 0.);
        pp.load("Delta", random_field_params.Delta, 1.);

        pp.load("calc_binned_power_spectrum", random_field_params.calc_binned_power_spectrum, 0);
        pp.load("bin_number", random_field_params.bin_number, random_field_params.N_readin/2); 
        pp.load("calc_higher_order_statistics", random_field_params.calc_higher_order_statistics, 0);
        pp.load("num_moments", random_field_params.num_orders, 0);
        pp.getarr("moments_to_print", random_field_params.orders, 0, random_field_params.num_orders);

        pp.load("read_from_stoiic", random_field_params.read_from_stoiic, 0);
        if(random_field_params.read_from_stoiic)
        {
            int num_modes;
            pp.load("n_k", num_modes, 0);
            pp.getarr("init_k", random_field_params.init_k, 0, num_modes);

            random_field_params.scalar_ps = amrex::Vector<amrex::Vector<amrex::Real>>(8, amrex::Vector<amrex::Real>(num_modes, 0.));
            pp.getarr("re_phi_k", random_field_params.scalar_ps[0], 0, num_modes);
            pp.getarr("im_phi_k", random_field_params.scalar_ps[1], 0, num_modes);
            pp.getarr("re_Pi_k", random_field_params.scalar_ps[2], 0, num_modes);
            pp.getarr("im_Pi_k", random_field_params.scalar_ps[3], 0, num_modes);
            pp.getarr("re_X_k", random_field_params.scalar_ps[4], 0, num_modes);
            pp.getarr("im_X_k", random_field_params.scalar_ps[5], 0, num_modes);
            pp.getarr("re_K_k", random_field_params.scalar_ps[6], 0, num_modes);
            pp.getarr("im_K_k", random_field_params.scalar_ps[7], 0, num_modes);

            random_field_params.tensor_ps = amrex::Vector<amrex::Vector<amrex::Real>>(4, amrex::Vector<amrex::Real>(num_modes, 0.));
            pp.getarr("re_h_k", random_field_params.tensor_ps[0], 0, num_modes);
            pp.getarr("im_h_k", random_field_params.tensor_ps[1], 0, num_modes);
            pp.getarr("re_dh_k", random_field_params.tensor_ps[2], 0, num_modes);
            pp.getarr("im_dh_k", random_field_params.tensor_ps[3], 0, num_modes);
        }
    }

    void check_params()
    {
        warn_parameter("scalar_mass", background_params.m,
                       background_params.m <
                           0.2 / coarsest_dx / dt_multiplier,
                       "oscillations of scalar field do not appear to be "
                       "resolved on coarsest level");

        warn_parameter("kstar", random_field_params.kstar,
                       random_field_params.kstar >= 0,
                       "cut-off frequency index must be positive");

        check_parameter("Delta", random_field_params.Delta,
                       (!random_field_params.calc_binned_power_spectrum
                        || random_field_params.Delta > 0),
                       "cut-off width must be positive and non-zero");

        check_parameter("orders", random_field_params.calc_higher_order_statistics,
                       (!random_field_params.calc_higher_order_statistics 
                        || !random_field_params.orders.empty()),
                       "moment orders must be provided");
    }

    // Initial data for matter and potential and BH
    double G_Newton;
    Potential::params_t potential_params;
    InitialBackgroundData::params_t background_params;
    InitialScalarData::params_t initial_params;
    RandomField::params_t random_field_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
