/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "BaseParameterChecker.hpp"
#include "GRParmParse.hpp"

class SimulationParameters
{
  public:
    SimulationParameters(GRParmParse &pp) = delete;

    static void check_params()
    {
       BaseParameterChecker::check_params();
       pp.load("num_points", num_points, 2);
    }

    int num_points;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
