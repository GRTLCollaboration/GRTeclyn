/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COORDINATES_HPP_
#define COORDINATES_HPP_

// Other includes
#include "DimensionDefinitions.hpp"
#include "simd.hpp"

#include <AMReX_IntVect.H>

#include <array>
#include <cmath>

template <class data_t> class Coordinates
{
  public:
    data_t x; // We vectorise over x so we must allow x to be a vector
    double y{};
    double z{};
    std::array<double, AMREX_SPACEDIM> m_center;

    AMREX_GPU_HOST_DEVICE
    Coordinates(amrex::IntVect integer_coords, double dx,
                std::array<double, AMREX_SPACEDIM> center = {0})
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
    compute_coord(double &out, int position, double dx,
                  double center_distance = 0)
    {
        out = (position + 0.5) * dx - center_distance;
    }

#if !defined(AMREX_USE_GPU)
    AMREX_FORCE_INLINE
    static typename std::enable_if_t<(simd_traits<double>::simd_len > 1), void>
    compute_coord(simd<double> &out, int position, double dx,
                  double center_distance = 0)
    {
        // NOLINTNEXTLINE
        double out_arr[simd_traits<double>::simd_len];
        for (int i = 0; i < simd_traits<double>::simd_len; ++i)
        {
            out_arr[i] = (position + i + 0.5) * dx - center_distance;
        }
        out = simd<double>::load(&out_arr[0]);
    }
#endif

    /// This function returns the radius subject to a floor for a given
    /// Coordinates object.
    AMREX_GPU_HOST_DEVICE [[nodiscard]] AMREX_FORCE_INLINE data_t
    get_radius() const
    {
        // Note that this is not currently dimension independent
        data_t r = sqrt(x * x + y * y + z * z);

        const double minimum_r = 1e-6;
        return simd_max(r, minimum_r);
    }

    /// This static function returns the radius subject to a floor
    /// for when no coordinates object exists.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static data_t
    get_radius(amrex::IntVect integer_coords, double dx,
               std::array<double, AMREX_SPACEDIM> center = {0})
    {
        data_t x;
        double y = NAN;
        double z = NAN;

        // Note that this is not currently dimension independent
        compute_coord(x, integer_coords[0], dx, center[0]);
        compute_coord(y, integer_coords[1], dx, center[1]);
        compute_coord(z, integer_coords[2], dx, center[2]);

        data_t r = std::sqrt(x * x + y * y + z * z);

        const double minimum_r = 1e-6;
        return simd_max(r, minimum_r);
    }
};

template <typename data_t>
ALWAYS_INLINE std::ostream &operator<<(std::ostream &a_os,
                                       const Coordinates<data_t> &in_coords)
{
    a_os << "(x,y,z) = (" << in_coords.x << "," << in_coords.y << ","
         << in_coords.z << ")"
         << " r = " << in_coords.get_radius();
    return a_os;
}
#endif /* COORDINATES_HPP_ */
