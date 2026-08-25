/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

#include "BaseParameterChecker.hpp"
#include "CCZ4RHS.hpp"
#include "MovingPunctureGauge.hpp"

#include <array>
#include <cmath>
#include <string>

class SimulationParameters
{
  public:
    SimulationParameters() = delete;

    static void check_params()
    {
        BaseParameterChecker::check_params();
        CCZ4_params_t::check_params();
        MovingPunctureGauge::params_t::check_params();

        GRParmParse scalar_field_pp("scalar_field");
        amrex::Real G_Newton{1.0};
        scalar_field_pp.queryAdd("G_Newton", G_Newton);
        if (G_Newton < 0.0)
        {
            scalar_field_pp.error("G_Newton", "must be >= 0.0");
        }

        amrex::Real scalar_mass{1.0};
        scalar_field_pp.queryAdd("scalar_mass", scalar_mass);
        if (scalar_mass < 0.0)
        {
            scalar_field_pp.error("scalar_mass", "must be >= 0.0");
        }

        GRParmParse geometry_pp("geometry");
        std::array<amrex::Real, AMREX_SPACEDIM> prob_extent{};
        geometry_pp.get("prob_extent", prob_extent);
        amrex::Real coarsest_dx{};
        geometry_pp.get("coarsest_dx", coarsest_dx);

        GRParmParse evolution_pp("evolution");
        amrex::Real dt_multiplier{};
        evolution_pp.get("dt_multiplier", dt_multiplier);
        if (scalar_mass >= 0.2 / coarsest_dx / dt_multiplier)
        {
            scalar_field_pp.warning(
                "scalar_mass", "oscillations of the scalar field may not be "
                               "resolved on the coarsest level");
        }

        GRParmParse line_pp("line_extraction");
        bool line_extraction_enabled{false};
        line_pp.queryAdd("enabled", line_extraction_enabled);
        std::string output_subpath{"extraction_data"};
        line_pp.queryAdd("output_subpath", output_subpath);

        if (line_extraction_enabled)
        {
            int num_points{128};
            line_pp.queryAdd("num_points", num_points);
            if (num_points < 2)
            {
                line_pp.error("num_points", "must be >= 2");
            }

            amrex::Real max_radius{amrex::Real(0.5) * prob_extent[0]};
            line_pp.queryAdd("max_radius", max_radius);
            if (max_radius <= 0.0)
            {
                line_pp.error("max_radius", "must be > 0.0");
            }

            std::array<amrex::Real, AMREX_SPACEDIM> center{};
            geometry_pp.get("center", center);
            const amrex::Real coordinate_offset =
                max_radius /
                std::sqrt(static_cast<amrex::Real>(AMREX_SPACEDIM));
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
            {
                if (center[dir] + coordinate_offset >= prob_extent[dir])
                {
                    line_pp.error(
                        "max_radius",
                        "places the line extraction outside the upper "
                        "boundary in at least one direction");
                }
            }
        }
    }
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
