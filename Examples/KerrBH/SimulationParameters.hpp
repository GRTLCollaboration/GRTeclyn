/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP
#define SIMULATIONPARAMETERS_HPP

// General includes
#include "BaseParameterChecker.hpp"

// Problem specific includes:
#include "CCZ4RHS.hpp"
#include "KerrBHInitialData.hpp"
#include "MovingPunctureGauge.hpp"
#include "PunctureTracker.hpp"

class SimulationParameters
{
  public:
    // NOLINTNEXTLINE(readability-identifier-length)
    SimulationParameters() = delete;

    static void check_params()
    {
        BaseParameterChecker::check_params();

        CCZ4_params_t::check_params();
        MovingPunctureGauge<FourthOrderDerivatives>::params_t::check_params();
        puncture_tracker_params_t::check_params();

        KerrBHInitialData::params_t::check_params();
    }
};

#endif /* SIMULATIONPARAMETERS_HPP */
