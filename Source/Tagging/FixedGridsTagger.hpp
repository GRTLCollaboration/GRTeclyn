/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef FIXEDGRIDSTAGGER_HPP_
#define FIXEDGRIDSTAGGER_HPP_

#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "GRParmParse.hpp"

// AMReX includes
#include <AMReX_TagBox.H>

using namespace amrex::literals;

class FixedGridsTagger
{
  protected:
    amrex::Real m_dx;
    int m_level;
    std::array<amrex::Real, AMREX_SPACEDIM> m_center{};
    std::array<amrex::Real, AMREX_SPACEDIM> m_domain_lengths{};

  public:
    static void check_params()
    {
        GRParmParse geometry_pp("geometry");
        std::array<amrex::Real, AMREX_SPACEDIM> center{};
        geometry_pp.get("center", center);
        std::array<amrex::Real, AMREX_SPACEDIM> domain_lengths{};
        geometry_pp.get("prob_extent", domain_lengths);

        for (int direction = 0; direction < AMREX_SPACEDIM; ++direction)
        {
            if (center[direction] < 0.0_rt ||
                center[direction] > domain_lengths[direction])
            {
                geometry_pp.error(
                    "center",
                    "must lie within the computational domain in every "
                    "direction when using FixedGridsTagger");
            }
        }
    }

    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    FixedGridsTagger(amrex::Real a_dx, int a_level)
        : m_dx(a_dx), m_level(a_level)
    {
        GRParmParse geometry_pp("geometry");
        geometry_pp.get("center", m_center);
        geometry_pp.get("prob_extent", m_domain_lengths);
    }

    AMREX_GPU_DEVICE void
    operator()(int ix, int iy, int iz,
               const amrex::Array4<amrex::TagBox::TagType> &tags) const
    // NOLINTEND(bugprone-easily-swappable-parameters)
    {
        // Refine half the simulated domain length on level 0, then halve the
        // refined length again on each subsequent level.
        const amrex::Real fraction_of_each_side_to_refine =
            std::pow(0.5_rt, m_level + 1.0_rt);
        const amrex::IntVect cell(AMREX_D_DECL(ix, iy, iz));

        for (int direction = 0; direction < AMREX_SPACEDIM; ++direction)
        {
            amrex::Real coordinate_from_center{};
            Coordinates::compute_coord(coordinate_from_center, cell[direction],
                                       m_dx, m_center[direction]);

            // clang-tidy incorrectly treats this initializer as missing.
            // NOLINTNEXTLINE(cppcoreguidelines-init-variables)
            amrex::Real distance_to_boundary =
                m_domain_lengths[direction] - m_center[direction];
            if (coordinate_from_center < 0.0_rt)
            {
                distance_to_boundary = m_center[direction];
            }

            const amrex::Real tagging_distance_from_center =
                fraction_of_each_side_to_refine * distance_to_boundary;
            if (std::abs(coordinate_from_center) >=
                tagging_distance_from_center)
            {
                return;
            }
        }

        tags(ix, iy, iz) = amrex::TagBox::SET;
    }
};

#endif /* FIXEDGRIDSTAGGER_HPP_ */
