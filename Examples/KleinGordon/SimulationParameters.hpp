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

    void read_params(GRParmParse &pp)
    {

        // These are parameters specfic to the Klein Gordon example

        pp.query("initial_amplitude", ampl);
        // What is their initial amplitude (if wave ICs )
        pp.query("initial_width", width);
        // What is the width of the Gaussian initial condition
        pp.query("scalar_mass",
                 scalar_mass); // What is the mass of the scalar particle?
        pp.query("wave_vector",
                 k_r);            // What is the wave number (if wave ICs)
        pp.query("alpha", alpha); // this is for Sine-Gordon ICs only
    }

    amrex::Real cfl = 0.2;
    amrex::Real ampl;
    amrex::Real width;
    int nfields             = 1;
    amrex::Real scalar_mass = 1.0;
    int ncomp               = 2;
    amrex::Real k_r         = 1.0;
    amrex::Real alpha       = 1.0;
    amrex::Real sigma       = 0.0;
};

#endif /* SIMULATIONPARAMETERS_HPP */
