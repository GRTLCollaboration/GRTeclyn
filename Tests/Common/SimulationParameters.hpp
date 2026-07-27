/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP
#define SIMULATIONPARAMETERS_HPP

// General includes
#include "GRParmParse.hpp"
#include "BaseParameterChecker.hpp"
#include "PunctureTracker.hpp"

class SimulationParameters : public BaseParameterChecker
{
  public:
    // NOLINTNEXTLINE(readability-identifier-length)
    SimulationParameters() : BaseParameterChecker()
    {
        check_params();
    }

    /// Read shared parameters
    // NOLINTNEXTLINE(readability-identifier-length)
    void check_params()
    {
        GRParmParse pp;
        bool puncture_tracking_enabled{false};
        pp.queryAdd("puncture_tracking.enabled", puncture_tracking_enabled);
        if (puncture_tracking_enabled)
        {
            puncture_tracker_params_t::check_params();
        }
        
        amrex::Real fake_bh1_mass = 0.5;
        amrex::Real fake_bh2_mass = 0.5;

        pp.queryAdd("fake_bh1_mass", fake_bh1_mass);
        pp.queryAdd("fake_bh2_mass", fake_bh2_mass);

        int num_points = 2;
        pp.queryAdd("num_points", num_points);

        bool verbosity = true;
        pp.queryAdd("verbosity", verbosity);
    }
};

#endif /* SIMULATIONPARAMETERS_HPP */
