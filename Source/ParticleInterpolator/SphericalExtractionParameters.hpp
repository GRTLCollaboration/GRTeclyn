/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SPHERICALEXTRACTIONPARAMETERS_HPP_
#define SPHERICALEXTRACTIONPARAMETERS_HPP_

// AMReX includes
#include <AMReX_GpuContainers.H>

// Parameters
#include "SurfaceExtractionParameters.hpp"

struct spherical_extraction_params_t : surface_extraction_params_t
{
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

    std::array<double, AMREX_SPACEDIM> center{}; //!< the center of the
                                                 //!< spherical shells
    std::array<double, AMREX_SPACEDIM> &extraction_center()
    {
        return this->center;
    }
    int num_modes{}; //!< the number of modes to extract
    std::vector<std::pair<int, int>> modes{}; //!< the modes to extract
                                              //!< l = first, m = second

    [[nodiscard]] const surface_extraction_params_t &
    get_surface_extraction_params() const
    {
        return *this;
    }

    static void check_params()
    {
        GRParmParse pp;
        int num_extraction_radii = 1;
        pp.queryAdd("num_extraction_radii", num_extraction_radii);
        if (num_extraction_radii <= 0)
        {
            pp.error("num_extraction_radii",
                     "must be bigger than 0 when activate_extraction = 1");
        }

        // Check for multiple extraction radii, otherwise load single
        // radius/level (for backwards compatibility).
        std::vector<int> extraction_levels(num_extraction_radii, 0);
        if (!pp.contains("extraction_levels"))
        {
            int extraction_level = 0;
            pp.queryAdd("extraction_level", extraction_level);
            std::fill(extraction_levels.begin(), extraction_levels.end(),
                      extraction_level);
            pp.addarr("extraction_levels", extraction_levels);
        }

        int min_extraction_level = *(std::min_element(extraction_levels.begin(),
                                                      extraction_levels.end()));
        pp.add("min_extraction_level", min_extraction_level);

        std::vector<double> extraction_radii(num_extraction_radii);
        if (!pp.contains("extraction_radii"))
        {
            double extraction_radius = 0.1;
            pp.queryAdd("extraction_radius", extraction_radius);
            std::fill(extraction_radii.begin(), extraction_radii.end(),
                      extraction_radius);
            pp.addarr("extraction_radii", extraction_radii);
        }

        int num_points_phi = 2;
        pp.queryAdd("num_points_phi", num_points_phi);

        int num_points_theta = 5;
        pp.queryAdd("num_points_theta", num_points_theta);
        if (num_points_theta % 2 == 0)
        {
            num_points_theta += 1;
            pp.warning("num_points_theta",
                       "incompatible with Simpson's rule so increased by 1");
            pp.add("num_points_theta", num_points_theta);
        }

        std::array<double, AMREX_SPACEDIM> amr_center{};
        pp.get("amr.center", amr_center);
        std::array<double, AMREX_SPACEDIM> extraction_center = amr_center;
        pp.queryAdd("extraction_center", extraction_center);

        // Used to check params

        std::array<int, AMREX_SPACEDIM> lo_boundary{};
        std::array<int, AMREX_SPACEDIM> hi_boundary{};
        std::array<double, AMREX_SPACEDIM> prob_extent{};

        pp.get("lo_boundary", lo_boundary);
        pp.get("hi_boundary", hi_boundary);
        pp.get("geometry.prob_extent", prob_extent);

        std::array<double, AMREX_SPACEDIM> reflective_domain_lo{};
        std::array<double, AMREX_SPACEDIM> reflective_domain_hi{};
        FOR (idir)
        {
            reflective_domain_lo[idir] =
                ((lo_boundary[idir] == BoundaryConditions::REFLECTIVE_BC)
                     ? -1.0
                     : 0.0) *
                prob_extent[idir];
            reflective_domain_hi[idir] =
                ((hi_boundary[idir] == BoundaryConditions::REFLECTIVE_BC)
                     ? 2.0
                     : 1.0) *
                prob_extent[idir];
        }

        FOR (idir)
        {
            // NOLINTNEXTLINE(cppcoreguidelines-init-variables)
            if (extraction_center[idir] < reflective_domain_lo[idir] ||
                extraction_center[idir] > reflective_domain_hi[idir])
            {
                pp.error("extraction_center",
                         "must be in the computational domain after applying "
                         "reflective symmetry");
            }

            for (int iradius = 0; iradius < num_extraction_radii; ++iradius)
            {
                // NOLINTNEXTLINE(cppcoreguidelines-init-variables)
                if (idir == 0)
                {
                    if (extraction_radii[iradius] < 0.0)
                    {
                        pp.error("extraction_radii", "must all be >= 0.0");
                    }
                }
                if (extraction_center[idir] - extraction_radii[iradius] <
                        reflective_domain_lo[idir] ||
                    extraction_center[idir] + extraction_radii[iradius] >
                        reflective_domain_hi[idir])
                {
                    pp.error(
                        "extraction_radii",
                        "extraction sphere must lie within the computational "
                        "domain after applying reflective symmetry");
                }
            }
        }

        int num_modes{};
        std::vector<int> extraction_modes_vect{};
        std::vector<std::pair<int, int>> modes{}; //!< the modes to extract
                                                  //!< l = first, m = second
        if (pp.contains("modes"))
        {
            pp.queryAdd("num_modes", num_modes);
            extraction_modes_vect.resize(static_cast<size_t>(2 * num_modes));
            pp.getarr("modes", extraction_modes_vect);
        }
        else
        {
            // by default extraction (l,m) = (2,0), (2,1) and (2,2)
            num_modes = 3;
            pp.add("num_modes", num_modes);
            extraction_modes_vect.resize(static_cast<size_t>(2 * num_modes));
            for (size_t i = 0; i < 3; ++i)
            {
                extraction_modes_vect[2 * i]     = 2;
                extraction_modes_vect[2 * i + 1] = static_cast<int>(i);
            }
            pp.addarr("modes", extraction_modes_vect);
        }

        modes.resize(num_modes);
        for (size_t i = 0; i < num_modes; ++i)
        {
            modes[i].first  = extraction_modes_vect[2 * i];
            modes[i].second = extraction_modes_vect[2 * i + 1];
        }

        for (int imode = 0; imode < num_modes; ++imode)
        {
            auto &mode = modes[imode];
            int l      = mode.first;
            int m      = mode.second;
            if ((l < 2) || (std::abs(m) > l))
            {
                pp.warning("modes",
                           "l must be >= 2 and m must satisfy -l <= m <= l");
            }
        }

        bool write_extraction = false;
        pp.queryAdd("write_extraction", write_extraction);

        std::string output_path = "";
        pp.get("output_path", output_path);

        std::string extraction_path = output_path+"/extraction_data/";
        std::string data_path = extraction_path;
        pp.add("extraction_path", extraction_path);
        pp.add("data_path", data_path);

        // default names to Weyl extraction
        std::string extraction_file_prefix = "Weyl4_extraction_";
        pp.queryAdd("extraction_file_prefix", extraction_file_prefix);
        std::string integral_file_prefix = "Weyl4_mode_";
        pp.queryAdd("integral_file_prefix", integral_file_prefix);
    };

    void fill_params()
    {
        GRParmParse pp;

        pp.get("num_extraction_radii", num_extraction_radii());

        // Check for multiple extraction radii, otherwise load single
        // radius/level (for backwards compatibility).
        std::vector<int> extraction_levels_stdvect(num_extraction_radii());
        pp.getarr("extraction_levels", extraction_levels_stdvect);

        extraction_levels.resize(num_extraction_radii());
        std::copy(extraction_levels_stdvect.begin(),
                  extraction_levels_stdvect.end(), extraction_levels.begin());

        std::vector<double> extraction_radii_stdvect(num_extraction_radii());
        pp.getarr("extraction_radii", extraction_radii_stdvect);

        extraction_radii().resize(num_extraction_radii());
        std::copy(extraction_radii_stdvect.begin(),
                  extraction_radii_stdvect.end(), extraction_radii().begin());

        pp.get("num_points_phi", num_points_phi());
        pp.get("num_points_theta", num_points_theta());

        pp.load("extraction_center", center);

        pp.get("num_modes", num_modes);
        std::vector<int> extraction_modes_vect{};
        pp.getarr("modes", extraction_modes_vect);

        modes.resize(num_modes);
        for (size_t i = 0; i < num_modes; ++i)
        {
            modes[i].first  = extraction_modes_vect[2 * i];
            modes[i].second = extraction_modes_vect[2 * i + 1];
        }

        pp.load("write_extraction", write_extraction);

        pp.get("data_path", data_path);
        pp.get("extraction_path", extraction_path);

        // default names to Weyl extraction
        pp.get("extraction_file_prefix", extraction_file_prefix);
        pp.get("integral_file_prefix", integral_file_prefix);
    };
};

#endif /* SPHERICALEXTRACTIONPARAMETERS_HPP_ */