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

#include <cmath>

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
        pp.load("activate_line_extraction", activate_line_extraction, false);
        if (activate_line_extraction)
        {
            pp.load("line_extraction_num_points", line_extraction_num_points,
                    128);
            pp.load("line_extraction_max_radius", line_extraction_max_radius,
                    0.5 * L);
        }
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
        if (activate_line_extraction)
        {
            check_parameter("line_extraction_num_points",
                            line_extraction_num_points,
                            line_extraction_num_points >= 2, "must be >= 2");
            check_parameter("line_extraction_max_radius",
                            line_extraction_max_radius,
                            line_extraction_max_radius > 0.0, "must be > 0.0");

            const amrex::Real coordinate_offset =
                line_extraction_max_radius /
                std::sqrt(static_cast<amrex::Real>(AMREX_SPACEDIM));
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
            {
                const amrex::Real end_coord = center[dir] + coordinate_offset;
                check_parameter(
                    "line extraction end coordinate", end_coord,
                    end_coord < reflective_domain_hi[dir],
                    "must lie inside the upper boundary in every direction");
            }
        }
    }

    amrex::Real G_Newton{1.0};
    bool activate_line_extraction{false};
    int line_extraction_num_points{128};
    amrex::Real line_extraction_max_radius{};
    OscillatonInitialData::params_t initial_params{};
    Potential::params_t potential_params{};
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
