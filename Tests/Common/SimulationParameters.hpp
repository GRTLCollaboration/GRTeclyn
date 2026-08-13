/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP
#define SIMULATIONPARAMETERS_HPP

// General includes
#include "BaseParameterChecker.hpp"
#include "BoostedBHInitialData.hpp"
#include "CCZ4RHS.hpp"
#include "GRParmParse.hpp"
#include "PunctureTagger.hpp"
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
        CCZ4_params_t::check_params();
        PunctureTagger<2>::check_params();
        puncture_tracker_params_t::check_params();

        // For the AHFinder test, which uses BinaryBH initial data
        BoostedBHInitialData::params_t::check_params(1);
        BoostedBHInitialData::params_t::check_params(2);

        GRParmParse test_pp("test");

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
