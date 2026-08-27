/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

#include "BaseParameterChecker.hpp"
#include "CCZ4RHS.hpp"
#include "FixedGridsTagger.hpp"
#include "LineExtractionParameters.hpp"
#include "MovingPunctureGauge.hpp"
#include "Potential.hpp"
#include "ScalarField.hpp"

class SimulationParameters
{
  public:
    SimulationParameters() = delete;

    static void check_params()
    {
        BaseParameterChecker::check_params();
        CCZ4_params_t::check_params();
        MovingPunctureGauge::params_t::check_params();
        FixedGridsTagger::check_params();

        ScalarField<Potential>::params_t::check_params();
        Potential::params_t::check_params();
        line_extraction_params_t::check_params("phi_line_extraction");
        line_extraction_params_t::check_params("rho_line_extraction");
    }
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
