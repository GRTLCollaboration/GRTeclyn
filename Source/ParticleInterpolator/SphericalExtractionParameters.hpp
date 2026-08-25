/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SPHERICALEXTRACTIONPARAMETERS_HPP_
#define SPHERICALEXTRACTIONPARAMETERS_HPP_

#include "BoundaryConditions.hpp"
#include "GRParmParse.hpp"
#include "SurfaceExtractionParameters.hpp"

#include <AMReX_GpuContainers.H>

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

struct spherical_extraction_params_t : surface_extraction_params_t
{
    bool enabled{};
    std::array<amrex::Real, AMREX_SPACEDIM> center{}; //!< center of the shells
    int num_modes{};                          //!< number of modes to extract
    std::vector<std::pair<int, int>> modes{}; //!< l = first, m = second

    explicit spherical_extraction_params_t(std::string a_param_scope)
        : m_param_scope(std::move(a_param_scope))
    {
    }

    int &num_extraction_radii() { return this->num_surfaces; }

    [[nodiscard]] const int &num_extraction_radii() const
    {
        return this->num_surfaces;
    }

    auto &extraction_radii() { return this->surface_param_values; }

    [[nodiscard]] const auto &extraction_radii() const
    {
        return this->surface_param_values;
    }

    int &num_points_theta() { return this->num_points_u; }

    [[nodiscard]] const int &num_points_theta() const
    {
        return this->num_points_u;
    }

    int &num_points_phi() { return this->num_points_v; }

    [[nodiscard]] const int &num_points_phi() const
    {
        return this->num_points_v;
    }

    [[nodiscard]] const surface_extraction_params_t &
    get_surface_extraction_params() const
    {
        return *this;
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

        GRParmParse geom_pp("geometry");
        std::array<amrex::Real, AMREX_SPACEDIM> grid_center{};
        geom_pp.get("center", grid_center);
        std::array<amrex::Real, AMREX_SPACEDIM> center = grid_center;
        extraction_pp.queryAdd("center", center);

        int num_radii = 1;
        extraction_pp.queryAdd("num_radii", num_radii);
        if (num_radii <= 0)
        {
            extraction_pp.error("num_radii", "must be greater than 0 when "
                                             "enabled = true");
        }

        std::vector<int> levels(num_radii, 0);
        if (extraction_pp.contains("levels"))
        {
            if (extraction_pp.countval("levels") != num_radii)
            {
                extraction_pp.error("levels", "must contain num_radii values");
            }
            extraction_pp.getarr("levels", levels, 0, num_radii);
        }
        else
        {
            extraction_pp.addarr("levels", levels);
        }

        const int min_level = *std::min_element(levels.begin(), levels.end());
        extraction_pp.add("min_level", min_level);

        std::vector<amrex::Real> radii(num_radii, 0.1);
        if (extraction_pp.contains("radii"))
        {
            if (extraction_pp.countval("radii") != num_radii)
            {
                extraction_pp.error("radii", "must contain num_radii values");
            }
            extraction_pp.getarr("radii", radii, 0, num_radii);
        }
        else
        {
            extraction_pp.addarr("radii", radii);
        }

        int num_points_phi = 2;
        extraction_pp.queryAdd("num_points_phi", num_points_phi);

        int num_points_theta = 5;
        extraction_pp.queryAdd("num_points_theta", num_points_theta);
        if (num_points_theta % 2 == 0)
        {
            num_points_theta += 1;
            extraction_pp.warning(
                "num_points_theta",
                "incompatible with Simpson's rule so increased by 1");
            extraction_pp.add("num_points_theta", num_points_theta);
        }

        std::array<int, AMREX_SPACEDIM> lo_condition{};
        std::array<int, AMREX_SPACEDIM> hi_condition{};
        std::array<amrex::Real, AMREX_SPACEDIM> prob_extent{};

        GRParmParse boundary_pp("boundary");
        lo_condition = BoundaryConditions::params_t::read_conditions(
            boundary_pp, "lo_condition");
        hi_condition = BoundaryConditions::params_t::read_conditions(
            boundary_pp, "hi_condition");
        geom_pp.get("prob_extent", prob_extent);

        std::array<amrex::Real, AMREX_SPACEDIM> reflective_domain_lo{};
        std::array<amrex::Real, AMREX_SPACEDIM> reflective_domain_hi{};
        FOR (idir)
        {
            reflective_domain_lo[idir] =
                ((lo_condition[idir] == BoundaryConditions::REFLECTIVE_BC)
                     ? -1.0
                     : 0.0) *
                prob_extent[idir];
            reflective_domain_hi[idir] =
                ((hi_condition[idir] == BoundaryConditions::REFLECTIVE_BC)
                     ? 2.0
                     : 1.0) *
                prob_extent[idir];
        }

        FOR (idir)
        {
            if (center[idir] < reflective_domain_lo[idir] ||
                center[idir] > reflective_domain_hi[idir])
            {
                extraction_pp.error("center",
                                    "must be in the computational domain after "
                                    "applying reflective symmetry");
            }

            for (int iradius = 0; iradius < num_radii; ++iradius)
            {
                if (idir == 0 && radii[iradius] < 0.0)
                {
                    extraction_pp.error("radii", "must all be >= 0.0");
                }
                if (center[idir] - radii[iradius] <
                        reflective_domain_lo[idir] ||
                    center[idir] + radii[iradius] > reflective_domain_hi[idir])
                {
                    extraction_pp.error("radii",
                                        "extraction sphere must lie within the "
                                        "computational domain after applying "
                                        "reflective symmetry");
                }
            }
        }

        int num_modes{};
        std::vector<int> modes_vector{};
        if (extraction_pp.contains("modes"))
        {
            extraction_pp.queryAdd("num_modes", num_modes);
            if (extraction_pp.countval("modes") != 2 * num_modes)
            {
                extraction_pp.error("modes", "must contain 2 * num_modes "
                                             "values");
            }
            modes_vector.resize(static_cast<size_t>(2) * num_modes);
            extraction_pp.getarr("modes", modes_vector);
        }
        else
        {
            // By default extract (l,m) = (2,0), (2,1) and (2,2).
            num_modes = 3;
            extraction_pp.add("num_modes", num_modes);
            modes_vector = {2, 0, 2, 1, 2, 2};
            extraction_pp.addarr("modes", modes_vector);
        }

        for (int imode = 0; imode < num_modes; ++imode)
        {
            const int l = modes_vector[2 * imode];
            const int m = modes_vector[2 * imode + 1];
            if ((l < 2) || (std::abs(m) > l))
            {
                extraction_pp.warning(
                    "modes", "l must be >= 2 and m must satisfy -l <= m <= l");
            }
        }

        bool write = false;
        extraction_pp.queryAdd("write", write);

        GRParmParse grteclyn_pp("grteclyn");
        std::string output_path{};
        grteclyn_pp.get("output_path", output_path);

        std::string path = output_path + "/extraction_data/";
        extraction_pp.queryAdd("path", path);

        std::string file_prefix = a_param_scope + "_";
        extraction_pp.queryAdd("file_prefix", file_prefix);
        std::string integral_file_prefix = a_param_scope + "_mode_";
        extraction_pp.queryAdd("integral_file_prefix", integral_file_prefix);
    }

    void fill_params()
    {
        GRParmParse extraction_pp(m_param_scope);

        extraction_pp.get("enabled", enabled);
        extraction_pp.get("center", center);
        if (!enabled)
        {
            return;
        }

        extraction_pp.get("num_radii", num_extraction_radii());

        std::vector<int> levels(num_extraction_radii());
        extraction_pp.getarr("levels", levels, 0, num_extraction_radii());
        extraction_levels.resize(num_extraction_radii());
        std::copy(levels.begin(), levels.end(), extraction_levels.begin());

        std::vector<amrex::Real> radii(num_extraction_radii());
        extraction_pp.getarr("radii", radii, 0, num_extraction_radii());
        extraction_radii().resize(num_extraction_radii());
        std::copy(radii.begin(), radii.end(), extraction_radii().begin());

        extraction_pp.get("num_points_phi", num_points_phi());
        extraction_pp.get("num_points_theta", num_points_theta());

        extraction_pp.get("num_modes", num_modes);
        std::vector<int> modes_vector(static_cast<size_t>(2 * num_modes));
        extraction_pp.getarr("modes", modes_vector);
        modes.resize(num_modes);
        for (int imode = 0; imode < num_modes; ++imode)
        {
            modes[imode].first  = modes_vector[2 * imode];
            modes[imode].second = modes_vector[2 * imode + 1];
        }

        extraction_pp.get("write", write_extraction);
        extraction_pp.get("path", data_path);
        extraction_path = data_path;
        extraction_pp.get("file_prefix", extraction_file_prefix);
        extraction_pp.get("integral_file_prefix", integral_file_prefix);
    }

  private:
    std::string m_param_scope;
};

#endif /* SPHERICALEXTRACTIONPARAMETERS_HPP_ */
