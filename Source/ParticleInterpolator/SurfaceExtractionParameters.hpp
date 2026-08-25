/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SURFACEEXTRACTIONPARAMETERS_HPP_
#define SURFACEEXTRACTIONPARAMETERS_HPP_

// AMReX includes
#include <AMReX_GpuContainers.H>

struct surface_extraction_params_t
{
    int num_surfaces{}; //!< number of surfaces over which to extraction
    amrex::Gpu::ManagedVector<amrex::Real>
        surface_param_values; //!< the values of the
                              //!< parameter that gives the required
                              //!< surfaces with SurfaceGeom geometry (e.g.
                              //!< radii for spherical shells)
    int num_points_u{};       //!< the number of points for the first parameter
                              //!< that parameterises each surface
    int num_points_v{};       //!< the number of points for the second parameter
                              //!< that parameterises each surfaces
    amrex::Gpu::ManagedVector<int>
        extraction_levels;   //!< the level on which to do the
                             //!< extraction for each surface
    bool write_extraction{}; //!< whether or not to write the extracted data

    std::string data_path, integral_file_prefix;
    std::string extraction_path, extraction_file_prefix;

    [[nodiscard]]
    int min_extraction_level() const
    {
        return *(std::min_element(extraction_levels.begin(),
                                  extraction_levels.end()));
    }

    void fill_params() {}
};

#endif /* SURFACEEXTRACTIONPARAMETERS_HPP_ */
