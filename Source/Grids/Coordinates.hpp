/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef COORDINATES_HPP_
#define COORDINATES_HPP_

// Other includes
#include "DimensionDefinitions.hpp"

#include <AMReX_IntVect.H>

#include <array>
#include <cmath>

class Coordinates
{
  public:
    amrex::Real x{};
    amrex::Real y{};
    amrex::Real z{};
    std::array<amrex::Real, AMREX_SPACEDIM> m_center;

    AMREX_GPU_HOST_DEVICE
    Coordinates(amrex::IntVect integer_coords, amrex::Real dx,
                std::array<amrex::Real, AMREX_SPACEDIM> center = {0.0})
        : m_center(center)
    {
        compute_coord(x, integer_coords[0], dx, center[0]);

// The below code allows for 2D Cartoon reduction:
#if DEFAULT_TENSOR_DIM == AMREX_SPACEDIM && AMREX_SPACEDIM == 3
        compute_coord(y, integer_coords[1], dx, center[1]);
        compute_coord(z, integer_coords[2], dx, center[2]);
#elif DEFAULT_TENSOR_DIM == AMREX_SPACEDIM + 1 && AMREX_SPACEDIM == 2
        y = 0;
        compute_coord(z, integer_coords[1], dx, center[1]);
#else
#ifdef AMREX_SPACEDIM
#error compute_coord has not got your dimension combination implemented.
#endif
#endif
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static void
    compute_coord(amrex::Real &out, int position, amrex::Real dx,
                  amrex::Real center_distance = 0.0)
    {
        out = (position + 0.5) * dx - center_distance;
    }

    /// This function returns the radius subject to a floor for a given
    /// Coordinates object.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
    get_radius() const
    {
        // Note that this is not currently dimension independent
        amrex::Real r = sqrt(x * x + y * y + z * z);
        return std::max(r, 1e-6);
    }

    /// This static function returns the radius subject to a floor
    /// for when no coordinates object exists.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static amrex::Real
    get_radius(amrex::IntVect integer_coords, amrex::Real dx,
               std::array<amrex::Real, AMREX_SPACEDIM> center = {0.0})
    {
        amrex::Real x = NAN;
        amrex::Real y = NAN;
        amrex::Real z = NAN;

        // Note that this is not currently dimension independent
        compute_coord(x, integer_coords[0], dx, center[0]);
        compute_coord(y, integer_coords[1], dx, center[1]);
        compute_coord(z, integer_coords[2], dx, center[2]);

        amrex::Real r = std::sqrt(x * x + y * y + z * z);
        return std::max(r, 1e-6);
    }
};

AMREX_FORCE_INLINE std::ostream &operator<<(std::ostream &a_os,
                                            const Coordinates &in_coords)
{
    a_os << "(x,y,z) = (" << in_coords.x << "," << in_coords.y << ","
         << in_coords.z << ")"
         << " r = " << in_coords.get_radius();
    return a_os;
}
#endif /* COORDINATES_HPP_ */
