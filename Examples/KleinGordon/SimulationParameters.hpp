/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP
#define SIMULATIONPARAMETERS_HPP

// General includes
#include "BaseParameterChecker.hpp"
#include "GRParmParse.hpp"

// Problem specific includes:
#include "ArrayTools.hpp"
#include "Wave.hpp"
#include "SineGordon.hpp"

class SimulationParameters : public BaseParameterChecker
{
  public:
    // NOLINTNEXTLINE(readability-identifier-length)
    SimulationParameters() : BaseParameterChecker()
    {
        // These parameters normally get read in inside SimulationParametersBase
        // but as this example doesn't use a lot of the other (CCZ4) parameters
        // this particular SimulationParameters class doesn't inherit from it.

        read_klein_gordon_params();
    }

    // NOLINTNEXTLINE(readability-identifier-length)
    void read_klein_gordon_params()
    {
        GRParmParse kg_pp("klein_gordon");

        // If the wave number isn't found in the params file
        // (so not wave ICs), look for the alpha parameter
        // (assume Sine-Gordon instead).

        std::string model{"Wave"};
        kg_pp.queryAdd("model", model);

        amrex::Real initial_time{0.0};
        kg_pp.queryAdd("initial_time", initial_time);

        if (model == "Wave")
        {
            Wave::params_t::check_params();
        }
        else if (model.find("SineGordon") == 0)
        {
            SineGordon::params_t::check_params();

        }
        else
        {
            amrex::Abort(
                "SimulationParameters: Klein Gordon model option not "
                "recognized. Choose from Wave, SineGordon1D or SineGordon3D.");
        }
    }
};

#endif /* SIMULATIONPARAMETERS_HPP */
