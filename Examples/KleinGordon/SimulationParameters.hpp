/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
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
        read_params(pp);
    }

    // NOLINTNEXTLINE(readability-identifier-length)
    void read_params(GRParmParse &pp)
    {

        // This parameter normally gets read in inside SimulationParametersBase
        // but this example doesn't use a lot of the other (CCZ4) parameters
        // so SimulationParameters doesn't inherit from it.
        pp.query("sigma", sigma);

        // If the wave number isn't found in the params file
        // (so not wave ICs), look for the alpha parameter
        // (assume Sine-Gordon instead).

        pp.query("model", model);

        if (model == "Wave")
        {

            pp.query("wave_vector", k_r);
            // Only wave example has the scalar mass as a parameter
            // SineGordon potential does not have a mass
            // associated with it.

            pp.query("scalar_mass",
                     scalar_mass); // What is the mass of the scalar particle?
        }
        else if (model.find("SineGordon") ==
                 0) // this is for Sine-Gordon ICs only
        {
            // These are parameters specfic to the Sine Gordon example
            pp.query("alpha", alpha);
        }
        else
        {
            amrex::Abort("Model option not recognized");
        }
    }

    amrex::Real scalar_mass{0.0};
    int ncomp{2};
    amrex::Real k_r{1.0};
    amrex::Real alpha{1.0};
    amrex::Real sigma{0.0};
    std::string model{"Wave"};
};

#endif /* SIMULATIONPARAMETERS_HPP */
