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
#include "BoostedBHInitialData.hpp"
#include "CCZ4RHS.hpp"
#include "PunctureTracker.hpp"
#ifdef USE_TWOPUNCTURES
#include "TwoPuncturesInitialData.hpp"
#endif

class SimulationParameters
{
  public:
    // NOLINTNEXTLINE(readability-identifier-length)
    SimulationParameters() = delete;

    static void check_params()
    {
        BaseParameterChecker::check_params();

        read_shared_params();

#ifndef USE_TWOPUNCTURES
        BoostedBHInitialData::params_t::check_params(1);
        BoostedBHInitialData::params_t::check_params(2);
#else
        TwoPuncturesInitialData::check_params();
#endif
    }

    /// Read shared parameters
    // NOLINTNEXTLINE(readability-identifier-length)
    static void read_shared_params()
    {
        GRParmParse pp;
        int formulation = CCZ4RHS<>::USE_CCZ4; // Whether to use BSSN or CCZ4
        pp.queryAdd("ccz4.formulation", formulation);

        if (formulation != CCZ4RHS<>::USE_CCZ4 &&
            formulation != CCZ4RHS<>::USE_BSSN)
        {
            pp.error("ccz4.formulation", "must be 0 or 1");
        }

        if (formulation == CCZ4RHS<>::USE_CCZ4)
        {
            CCZ4_params_t::check_params();
        }
        else if (formulation == CCZ4RHS<>::USE_BSSN)
        {
            if (pp.contains("ccz4.kappa1") || pp.contains("ccz4.kappa2") ||
                pp.contains("ccz4.kappa3"))
            {
                pp.warning("kappa1/2/3",
                           "should not be provided with BSSN formulation, "
                           "setting them all to zero");
            }
            pp.add("ccz4.kappa1", 0.0);
            pp.add("ccz4.kappa2", 0.0);
            pp.add("ccz4.kappa3", 0.0);
        }

        // Do we want puncture tracking?
        bool puncture_tracking_enabled{false};
        pp.queryAdd("puncture_tracking.enabled", puncture_tracking_enabled);
        if (puncture_tracking_enabled)
        {
            puncture_tracker_params_t::check_params();
        }
    }
};

#endif /* SIMULATIONPARAMETERS_HPP */
