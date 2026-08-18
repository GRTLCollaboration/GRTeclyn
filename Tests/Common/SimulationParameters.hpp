/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP
#define SIMULATIONPARAMETERS_HPP

// General includes
#include "BaseParameterChecker.hpp"
#include "CCZ4RHS.hpp"
#include "GRParmParse.hpp"
#include "PunctureTracker.hpp"
#include "SphericalExtractionParameters.hpp"

class SimulationParameters
{
  public:
    // NOLINTNEXTLINE(readability-identifier-length)
    SimulationParameters() = delete;

    /// Read shared parameters
    // NOLINTNEXTLINE(readability-identifier-length)
    static void check_params()
    {
        BaseParameterChecker::check_params();

        spherical_extraction_params_t::check_params("test_extraction_lo");
        spherical_extraction_params_t::check_params("test_extraction_hi");

        GRParmParse pp;
        GRParmParse test_pp("test");
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

        bool puncture_tracking_enabled{false};
        pp.queryAdd("puncture_tracking.enabled", puncture_tracking_enabled);
        if (puncture_tracking_enabled)
        {
            puncture_tracker_params_t::check_params();
        }

        amrex::Real fake_bh1_mass = 0.5;
        amrex::Real fake_bh2_mass = 0.5;

        test_pp.queryAdd("fake_bh1_mass", fake_bh1_mass);
        test_pp.queryAdd("fake_bh2_mass", fake_bh2_mass);

        int num_points = 2;
        test_pp.queryAdd("num_points", num_points);

        int es = 0;
        int el = 2;
        int em = 0;
        test_pp.queryAdd("es", es);
        test_pp.queryAdd("el", el);
        test_pp.queryAdd("em", em);
    }
};

#endif /* SIMULATIONPARAMETERS_HPP */
