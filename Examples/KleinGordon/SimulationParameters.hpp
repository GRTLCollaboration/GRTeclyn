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

        /* These used to be in KleinGordonLevel ... */

        // ParmParse pp("wave");
        //    pp.query("v", verbose); // Could use this to control verbosity
        //    during the run
        pp.query("nfields", nfields);
        pp.getarr("initial_amplitude", ampl, 0, nfields);
        pp.getarr("initial_width", width, 0, nfields);
        pp.query("scalar_mass", scalar_mass);
        pp.query("wave_vector", k_r);
        pp.query("alpha", alpha);

        ncomp = 2 * nfields;
    }

    amrex::Real cfl = 0.2;
    amrex::Vector<float> ampl;
    amrex::Vector<float> width;
    int nfields             = 1;
    amrex::Real scalar_mass = 1.0;
    int ncomp               = nfields * 2;
    amrex::Real k_r         = 1.0;
    amrex::Real alpha       = 1.0;
};

#endif /* SIMULATIONPARAMETERS_HPP */
