/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef PUNCTURETAGGER_HPP_
#define PUNCTURETAGGER_HPP_

#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"

#include <AMReX_Array4.H>
#include <AMReX_TagBox.H>

//! This class tags cells near the punctures so that the BH apparent horizons
//! are covered
template <unsigned int num_punctures> class PunctureTagger
{
  protected:
    amrex::Real m_dx;
    int m_level;
    int m_max_level;
    static constexpr unsigned int num_puncture_coords =
        AMREX_SPACEDIM * num_punctures;
    std::array<amrex::Real, num_punctures> m_puncture_masses;
    std::array<amrex::Real, num_puncture_coords> m_puncture_coords;

  public:

    // The constructor
    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    PunctureTagger(
        const amrex::Real a_dx, const int a_level, const int a_max_level,
        const std::array<amrex::Real, num_puncture_coords> &a_puncture_coords,
        const std::array<amrex::Real, num_punctures> &a_puncture_masses)
        // NOLINTEND(bugprone-easily-swappable-parameters)
        : m_dx(a_dx), m_level(a_level), m_max_level(a_max_level),
          m_puncture_masses(a_puncture_masses),
          m_puncture_coords(a_puncture_coords) {};

    AMREX_GPU_DEVICE void
    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    operator()(int i, int j, int k,
               const amrex::Array4<amrex::TagBox::TagType> &tags) const
    // NOLINTEND(bugprone-easily-swappable-parameters)
    {
        // ensure that the horizons of the punctures are covered
        // by the max level - for this we need
        // only check the puncture locations on the top 2 levels
        // which regrid (ie, max_level - 1 to max_level - 2)
        // (just the top level would be ok, but doing two ensures
        // the top levels are well spaced)

        // we want each level to be double the innermost one in size
        const int exponent  = std::min(m_max_level - m_level - 1, 1);
        const double factor = std::pow(2.0, exponent);

        amrex::IntVect current_cell(AMREX_D_DECL(i, j, k));
        // loop over puncture masses
        for (int ipuncture = 0; ipuncture < num_punctures; ++ipuncture)
        {
            std::array<amrex::Real, AMREX_SPACEDIM> current_puncture_coords = {
                AMREX_D_DECL(
                    m_puncture_coords[ipuncture * AMREX_SPACEDIM + 0],
                    m_puncture_coords[ipuncture * AMREX_SPACEDIM + 1],
                    m_puncture_coords[ipuncture * AMREX_SPACEDIM + 2])};

            const Coordinates<amrex::Real> coords(current_cell, m_dx,
                                                  current_puncture_coords);
            const amrex::Real r = coords.get_radius();
            // decide whether to tag based on distance to horizon
            // plus a fudge factor of 1.5
            amrex::Real fudge_factor = 1.5;
            if (r < fudge_factor * factor * m_puncture_masses[ipuncture])
            {
                tags(current_cell) = amrex::TagBox::SET;
            }
        }
    }
};

#endif /* PUNCTURETAGGER_HPP_ */
