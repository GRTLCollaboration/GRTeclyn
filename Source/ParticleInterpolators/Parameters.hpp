#ifndef PARAMETERS_HPP_
#define PARAMETERS_HPP_

#include <AMReX.H>   // for amrex::Real and AMREX_SPACEDIM
#include <algorithm> // for std::min_element
#include <array>
#include <string>
#include <vector>

struct Spherical_params_t
{
    int num_extraction_radii; //!< number of surfaces over which to extract
    std::array<amrex::Real, AMREX_SPACEDIM> extraction_center{};
    amrex::Real extraction_radii; //!< parameter for surfaces (e.g., radii)
    int num_points_theta;
    int num_points_phi;
    // std::vector<int> extraction_levels; //!< extraction level for each
    // surface
    bool write_extraction    = true;
    std::string chi_filename = "chi_spherical";

    // int min_extraction_level() const
    // {
    //     return *(std::min_element(extraction_levels.begin(),
    //     extraction_levels.end()));
    // }
};

#endif // PARAMETERS_HPP_
