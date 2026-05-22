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
    int num_modes{};                        //!< the number of modes to extract
    std::vector<std::pair<int, int>> modes; //!< the modes to extract
                                            //!< l = first, m = second

    [[nodiscard]] const surface_extraction_params_t &
    get_surface_extraction_params() const
    {
        return *this;
    }
};

#endif /* SPHERICALEXTRACTIONPARAMETERS_HPP_ */