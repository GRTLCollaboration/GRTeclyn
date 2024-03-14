/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "AMReXParameters.hpp"
#include "GRParmParse.hpp"

class SimulationParameters : public AMReXParameters
{
  public:
    SimulationParameters(GRParmParse &pp) : AMReXParameters(pp)
    {
        pp.load("num_points", num_points, 2);
    }

    int num_points;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
