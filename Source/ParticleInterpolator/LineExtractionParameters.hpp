/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef LINEEXTRACTIONPARAMETERS_HPP_
#define LINEEXTRACTIONPARAMETERS_HPP_

#include "BoundaryConditions.hpp"
#include "GRParmParse.hpp"

#include <array>
#include <cmath>
#include <string>
#include <utility>

struct line_extraction_params_t
{
    bool enabled{};
    int num_points{};
    std::array<amrex::ParticleReal, AMREX_SPACEDIM> start_coords{};
    std::array<amrex::ParticleReal, AMREX_SPACEDIM> end_coords{};
    std::string data_path;
    std::string file_prefix;

    explicit line_extraction_params_t(std::string a_param_scope)
        : m_param_scope(std::move(a_param_scope))
    {
    }

    static void check_params(const std::string &a_param_scope)
    {
        GRParmParse extraction_pp(a_param_scope);

        bool enabled = false;
        extraction_pp.queryAdd("enabled", enabled);
        if (!enabled)
        {
            return;
        }

        int num_points = 128;
        extraction_pp.queryAdd("num_points", num_points);
        if (num_points < 2)
        {
            extraction_pp.error("num_points", "must be >= 2");
        }

        GRParmParse geometry_pp("geometry");
        std::array<amrex::Real, AMREX_SPACEDIM> start_coords{};
        geometry_pp.get("center", start_coords);
        extraction_pp.queryAdd("start", start_coords);

        std::array<amrex::Real, AMREX_SPACEDIM> direction{};
        direction.fill(1.0);
        extraction_pp.queryAdd("direction", direction);
        amrex::Real direction_norm_squared = 0.0;
        for (const auto direction_component : direction)
        {
            direction_norm_squared += direction_component * direction_component;
        }
        if (direction_norm_squared == 0.0)
        {
            extraction_pp.error("direction", "must be nonzero");
        }

        std::array<amrex::Real, AMREX_SPACEDIM> prob_extent{};
        geometry_pp.get("prob_extent", prob_extent);
        amrex::Real max_radius = 0.5 * prob_extent[0];
        extraction_pp.queryAdd("max_radius", max_radius);
        if (max_radius <= 0.0)
        {
            extraction_pp.error("max_radius", "must be > 0.0");
        }

        const amrex::Real direction_norm = std::sqrt(direction_norm_squared);
        std::array<amrex::Real, AMREX_SPACEDIM> end_coords{};
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
        {
            end_coords[dir] = start_coords[dir] +
                              max_radius * direction[dir] / direction_norm;
        }

        GRParmParse boundary_pp("boundary");
        const auto lo_condition = BoundaryConditions::params_t::read_conditions(
            boundary_pp, "lo_condition");
        const auto hi_condition = BoundaryConditions::params_t::read_conditions(
            boundary_pp, "hi_condition");
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
        {
            const amrex::Real domain_lo =
                (lo_condition[dir] == BoundaryConditions::REFLECTIVE_BC)
                    ? -prob_extent[dir]
                    : 0.0;
            const amrex::Real domain_hi =
                (hi_condition[dir] == BoundaryConditions::REFLECTIVE_BC)
                    ? 2.0 * prob_extent[dir]
                    : prob_extent[dir];
            if (start_coords[dir] < domain_lo ||
                start_coords[dir] > domain_hi || end_coords[dir] < domain_lo ||
                end_coords[dir] > domain_hi)
            {
                extraction_pp.error(
                    "max_radius",
                    "places the line outside the computational domain after "
                    "applying reflective symmetry");
            }
        }

        GRParmParse grteclyn_pp("grteclyn");
        std::string output_path{};
        grteclyn_pp.get("output_path", output_path);
        std::string path = output_path + "/extraction_data/";
        extraction_pp.queryAdd("path", path);

        std::string file_prefix = a_param_scope + "_";
        extraction_pp.queryAdd("file_prefix", file_prefix);
    }

    void fill_params()
    {
        GRParmParse extraction_pp(m_param_scope);
        extraction_pp.get("enabled", enabled);
        if (!enabled)
        {
            return;
        }

        extraction_pp.get("num_points", num_points);
        extraction_pp.get("start", start_coords);

        std::array<amrex::ParticleReal, AMREX_SPACEDIM> direction{};
        extraction_pp.get("direction", direction);
        amrex::ParticleReal max_radius{};
        extraction_pp.get("max_radius", max_radius);
        amrex::ParticleReal direction_norm_squared = 0.0;
        for (const auto direction_component : direction)
        {
            direction_norm_squared += direction_component * direction_component;
        }
        const amrex::ParticleReal direction_norm =
            std::sqrt(direction_norm_squared);
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
        {
            end_coords[dir] = start_coords[dir] +
                              max_radius * direction[dir] / direction_norm;
        }

        extraction_pp.get("path", data_path);
        extraction_pp.get("file_prefix", file_prefix);
    }

  private:
    std::string m_param_scope;
};

#endif /* LINEEXTRACTIONPARAMETERS_HPP_ */
