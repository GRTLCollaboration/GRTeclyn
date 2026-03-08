/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP
#define SIMULATIONPARAMETERS_HPP

// General includes
#include "AMReXParameters.hpp"
#include "GRParmParse.hpp"

// Problem specific includes:
#include "ArrayTools.hpp"

class SimulationParameters : public AMReXParameters
{
  public:
    // NOLINTNEXTLINE(readability-identifier-length)
    SimulationParameters(GRParmParse &pp) : AMReXParameters(pp)
    {
        // These parameters normally get read in inside SimulationParametersBase
        // but as this example doesn't use a lot of the other (CCZ4) parameters
        // this particular SimulationParameters class doesn't inherit from it.

        pp.queryAdd("sigma", sigma);

        pp.queryAdd("nan_check", nan_check);

        read_klein_gordon_params(pp);
    }

    // NOLINTNEXTLINE(readability-identifier-length)
    void read_klein_gordon_params(GRParmParse &pp)
    {
        // If the wave number isn't found in the params file
        // (so not wave ICs), look for the alpha parameter
        // (assume Sine-Gordon instead).

        pp.queryAdd("klein_gordon.model", model);
        pp.queryAdd("klein_gordon.initial_time", t0);

        if (model == "Wave")
        {

            pp.queryAdd("klein_gordon.wave_vector", k_r);
            // Only wave example has the scalar mass as a parameter
            // SineGordon potential does not have a mass
            // associated with it.

            pp.queryAdd(
                "klein_gordon.scalar_mass",
                scalar_mass); // What is the mass of the scalar particle?
        }
        else if (model.find("SineGordon") ==
                 0) // this is for Sine-Gordon ICs only
        {
            // These are parameters specfic to the Sine Gordon example
            pp.queryAdd("klein_gordon.alpha", alpha);
        }
        else
        {
            amrex::Abort(
                "SimulationParameters: Klein Gordon model option not "
                "recognized. Choose from Wave, SineGordon1D or SineGordon3D.");
        }
    }

    static const int ncomp{2};

    amrex::Real alpha{1.0};
    amrex::Real k_r{1.0};
    amrex::Real scalar_mass{0.0};
    amrex::Real sigma{0.0};
    amrex::Real t0{0.0};

    std::string model{"Wave"};
    bool nan_check{true};
};

#endif /* SIMULATIONPARAMETERS_HPP */
