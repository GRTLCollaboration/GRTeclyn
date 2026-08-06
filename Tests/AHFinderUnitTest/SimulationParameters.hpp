/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP
#define SIMULATIONPARAMETERS_HPP

// General includes
#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"

// Problem specific includes: reuses the BinaryBH example's initial data
#include "BoostedBHInitialData.hpp"

class SimulationParameters : public SimulationParametersBase
{
  public:
    // NOLINTNEXTLINE(readability-identifier-length)
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        read_bh_params(pp);
    }

    /// Read BH parameters (same param names as the BinaryBH example)
    // NOLINTNEXTLINE(readability-identifier-length)
    void read_bh_params(GRParmParse &pp)
    {
        pp.load("massA", bh1_params.mass);
        pp.load("momentumA", bh1_params.momentum, {0.0, 0.0, 0.0});
        pp.load("massB", bh2_params.mass);
        pp.load("momentumB", bh2_params.momentum, {0.0, 0.0, 0.0});

        // Get the centers of the BHs either explicitly or as an offset from
        // the domain center (not both, or they will be offset from the
        // center provided)
        std::array<double, AMREX_SPACEDIM> centerA{};
        std::array<double, AMREX_SPACEDIM> centerB{};
        std::array<double, AMREX_SPACEDIM> offsetA{};
        std::array<double, AMREX_SPACEDIM> offsetB{};
        pp.load("centerA", centerA, center);
        pp.load("centerB", centerB, center);
        pp.load("offsetA", offsetA, {0.0, 0.0, 0.0});
        pp.load("offsetB", offsetB, {0.0, 0.0, 0.0});
        FOR (idir)
        {
            bh1_params.center[idir] = centerA[idir] + offsetA[idir];
            bh2_params.center[idir] = centerB[idir] + offsetB[idir];
        }
    }

    // Collection of parameters necessary for the BinaryBH initial data
    BoostedBHInitialData::params_t bh1_params{};
    BoostedBHInitialData::params_t bh2_params{};
};

#endif /* SIMULATIONPARAMETERS_HPP */
