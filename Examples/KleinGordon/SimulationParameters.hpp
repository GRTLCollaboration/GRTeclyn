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

        // These are parameters specfic to the Klein Gordon example

        pp.query("scalar_mass",
                 scalar_mass); // What is the mass of the scalar particle?

        // If the wave number isn't found in the params file
        // (so not wave ICs), look for the alpha parameter
        // (assume Sine-Gordon instead).
        if (!pp.query("wave_vector", k_r))
        {
            pp.query("alpha", alpha); // this is for Sine-Gordon ICs only
            model = "SineGordon";
        }
        else
        {
            model = "Wave";
        }
        pp.add("model", model);
    }

    amrex::Real cfl{0.2};
    amrex::Real scalar_mass{1.0};
    int ncomp{2};
    amrex::Real k_r{1.0};
    amrex::Real alpha{1.0};
    amrex::Real sigma{0.0};
    std::string model{"Wave"};
};

#endif /* SIMULATIONPARAMETERS_HPP */
