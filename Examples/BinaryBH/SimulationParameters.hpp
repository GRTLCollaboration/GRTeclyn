/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP
#define SIMULATIONPARAMETERS_HPP

// General includes
#include "BaseParameterChecker.hpp"

// Problem specific includes:
#include "ArrayTools.hpp"
#include "BoostedBHInitialData.hpp"
#include "CCZ4RHS.hpp"
#include "PunctureTracker.hpp"
#include "SphericalExtractionParameters.hpp"
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

        CCZ4_params_t::check_params();
        puncture_tracker_params_t::check_params();

#ifndef USE_TWOPUNCTURES
        BoostedBHInitialData::params_t::check_params(1);
        BoostedBHInitialData::params_t::check_params(2);
#else
        TwoPuncturesInitialData::check_params();
#endif

        spherical_extraction_params_t::check_params("weyl_extraction");
    }
};

#endif /* SIMULATIONPARAMETERS_HPP */
