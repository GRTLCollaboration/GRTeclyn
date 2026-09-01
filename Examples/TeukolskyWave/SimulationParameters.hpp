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
#include <utility>
#include <string>

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

        // Check the Teukolsky wave specific parameters
        GRParmParse tw_pp("teukolsky_wave");

        int magnetic{0};
        std::string parity{"even"};
        tw_pp.queryAdd("magnetic", magnetic);
        tw_pp.queryAdd("parity", parity);
        const std::set<std::pair<std::string, int>> implemented_combinations = {
            {"even", 0},
            {"even", 2},
            {"odd",  2}
        };
        if (implemented_combinations.find({parity, magnetic}) ==
            implemented_combinations.end())
        {
            tw_pp.error("Combination of magnetic/parity not implemented.  Must "
                        "be (even, 0), (even, 2) or (odd, 2).");
        }
        else
        {
            TeukolskyWaveInitialData::params_t::check_params();
        }
    }
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
