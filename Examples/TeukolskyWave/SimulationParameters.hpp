/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

#include "BaseParameterChecker.hpp"
#include "CCZ4RHS.hpp"
#include "FixedGridsTagger.hpp"
#include "IntegratedMovingPunctureGauge.hpp"
#include "SphericalExtractionParameters.hpp"
#include "TeukolskyWaveInitialData.hpp"

#include <set>
#include <string>
#include <utility>

class SimulationParameters
{
  public:
    SimulationParameters() = delete;

    static void check_params()
    {
        BaseParameterChecker::check_params();
        CCZ4_params_t::check_params();
        IntegratedMovingPunctureGauge<
            FourthOrderDerivatives>::params_t::check_params();
        FixedGridsTagger::check_params();

        spherical_extraction_params_t::check_params("weyl_extraction");

	TeukolskyWaveInitialData::params_t::check_params();
    }
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
