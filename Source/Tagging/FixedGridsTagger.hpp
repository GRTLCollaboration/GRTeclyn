/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef FIXEDGRIDSTAGGER_HPP_
#define FIXEDGRIDSTAGGER_HPP_

#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"

// AMReX includes
#include <AMReX_TagBox.H>

class FixedGridsTagger
{
  protected:
    double m_dx;
    double m_L;
    int m_level;
    std::array<double, AMREX_SPACEDIM> m_center;

  public:
    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    FixedGridsTagger(const double dx, const int a_level, const double a_L,
                     const std::array<double, AMREX_SPACEDIM> a_center)
        : m_dx(dx), m_L(a_L), m_level(a_level), m_center(a_center) {};
    // NOLINTEND(bugprone-easily-swappable-parameters)
    AMREX_GPU_DEVICE void
    operator()(int ix, int iy, int iz,
               const amrex::Array4<amrex::TagBox::TagType> &tags) const
    {
        // make sure the inner part is regridded around the horizon
        // take L as the length of full grid, so tag inner 1/2
        // of it, which means inner \pm L/4
        double ratio = pow(2.0, -(m_level + 2.0));

        amrex::IntVect cell(AMREX_D_DECL(ix, iy, iz));

        const Coordinates coords(cell, m_dx, m_center);
        const amrex::Real max_abs_xy =
            std::max(std::abs(coords.x), std::abs(coords.y));
        const amrex::Real max_abs_xyz =
            std::max(max_abs_xy, std::abs(coords.z));

        if (max_abs_xyz < m_L * ratio)
        {
            tags(ix, iy, iz) = amrex::TagBox::SET;
        }
    }
};

#endif /* FIXEDGRIDSTAGGER_HPP_ */
