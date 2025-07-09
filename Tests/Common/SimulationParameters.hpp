/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP
#define SIMULATIONPARAMETERS_HPP

// General includes
#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"

class SimulationParameters : public SimulationParametersBase
{
  public:
    // NOLINTNEXTLINE(readability-identifier-length)
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        read_params(pp);
    }

    /// Read shared parameters
    // NOLINTNEXTLINE(readability-identifier-length)
    void read_params(GRParmParse &pp)
    {
        int max_level = -1;
        pp.get("amr.max_level", max_level);
        // Do we want puncture tracking and constraint norm calculation?
        pp.load("puncture_tracking.enabled", puncture_tracking_enabled, true);
        pp.load("puncture_tracking.level", puncture_tracking_level, max_level);
        pp.load("puncture_tracking.initial_coords",
                puncture_tracking_initial_coords,
                {center[0], center[1] - 1.0, center[2], center[0],
                 center[1] + 1.0, center[2]});

        pp.load("fake_bh1_mass", fake_bh1_mass, 0.5);
        pp.load("fake_bh2_mass", fake_bh2_mass, 0.5);
    }

    bool puncture_tracking_enabled{};
    int puncture_tracking_level{};
    std::array<amrex::Real, AMREX_SPACEDIM * 2UL>
        puncture_tracking_initial_coords{};

    amrex::Real fake_bh1_mass{};
    amrex::Real fake_bh2_mass{};
};

#endif /* SIMULATIONPARAMETERS_HPP */
