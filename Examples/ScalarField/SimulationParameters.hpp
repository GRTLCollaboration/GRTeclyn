/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"

#include "OscillatonInitialData.hpp"
#include "Potential.hpp"

class SimulationParameters : public SimulationParametersBase
{
  public:
    // NOLINTNEXTLINE(readability-identifier-length)
    explicit SimulationParameters(GRParmParse &pp)
        : SimulationParametersBase(pp)
    {
        read_scalar_field_params(pp);
        check_scalar_field_params();
    }

    // NOLINTNEXTLINE(readability-identifier-length)
    void read_scalar_field_params(GRParmParse &pp)
    {
        initial_params.center = center;
        pp.load("G_Newton", G_Newton, 1.0);
        pp.load("scalar_mass", potential_params.scalar_mass, 1.0);
    }

    void check_scalar_field_params()
    {
        check_parameter("G_Newton", G_Newton, G_Newton >= 0.0,
                        "must be >= 0.0");
        check_parameter("scalar_mass", potential_params.scalar_mass,
                        potential_params.scalar_mass >= 0.0, "must be >= 0.0");
        warn_parameter("scalar_mass", potential_params.scalar_mass,
                       potential_params.scalar_mass <
                           0.2 / coarsest_dx / dt_multiplier,
                       "oscillations of the scalar field may not be resolved "
                       "on the coarsest level");
    }

    amrex::Real G_Newton{1.0};
    OscillatonInitialData::params_t initial_params{};
    Potential::params_t potential_params{};
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
