/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef EXTRACTIONTAGGER_HPP_
#define EXTRACTIONTAGGER_HPP_

#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "GRParmParse.hpp"
#include "SphericalExtractionParameters.hpp"
#include "Tensor.hpp"

#include <AMReX_Array4.H>
#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>
#include <AMReX_TagBox.H>

class ExtractionTagger
{
  protected:
    amrex::Real m_dx;
    int m_num_extraction_radii{};
    const amrex::Real *m_extraction_radii_ptr{nullptr};
    const int *m_extraction_levels_ptr{nullptr};
    std::array<amrex::Real, AMREX_SPACEDIM> m_center;
    int m_level;
    amrex::Real m_level_separation{1.5};

  public:
    static void check_params()
    {
        GRParmParse extraction_tagging_pp("extraction_tagging");
        amrex::Real level_separation{1.5};
        extraction_tagging_pp.queryAdd("level_separation", level_separation);

        if (level_separation < 1.0)
        {
            extraction_tagging_pp.warning(
                "level_separation",
                "levels may be too close together, which results in boundary "
                "error reflecting; either increase this value or set n_proper "
                "to be larger");
        }
        if (level_separation > 2.0)
        {
            extraction_tagging_pp.warning(
                "level_separation",
                "levels are more than doubling around punctures, which may "
                "result in too much refinement");
        }
    }

    // a_params must outlive the tagger and any GPU kernel which captures it.
    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
    ExtractionTagger(const amrex::Real dx, const int a_level,
                     const spherical_extraction_params_t &a_params)
        : m_dx(dx), m_level(a_level)
    {
        m_center = a_params.center;
        if (a_params.enabled)
        {
            m_num_extraction_radii  = a_params.num_extraction_radii();
            m_extraction_radii_ptr  = a_params.extraction_radii().data();
            m_extraction_levels_ptr = a_params.extraction_levels.data();
        }
        else
        {
            // Avoids conditionals in the kernel by setting num to 0
            m_num_extraction_radii = 0;
        }
        GRParmParse extraction_tagging_pp("extraction_tagging");
        extraction_tagging_pp.get("level_separation", m_level_separation);
    }

    AMREX_GPU_DEVICE void
    operator()(int i, int j, int k,
               const amrex::Array4<amrex::TagBox::TagType> &tags) const
    {
        // Enforce a given resolution at the extraction radii
        amrex::IntVect cell(AMREX_D_DECL(i, j, k));
        for (int iradius = 0; iradius < m_num_extraction_radii; ++iradius)
        {
            // Regrid if current level is not the required refinement level
            if (m_level < m_extraction_levels_ptr[iradius])
            {
                const Coordinates coords(cell, m_dx, m_center);
                const amrex::Real r = coords.get_radius();

                // Keep the levels spaced out
                const int exponent =
                    m_extraction_levels_ptr[iradius] - m_level - 1;
                const amrex::Real factor =
                    std::pow(m_level_separation, exponent);
                // Add a 20% buffer to extraction zone so not too near boundary
                if (r < 1.2 * factor * m_extraction_radii_ptr[iradius])
                {
                    tags(i, j, k) = amrex::TagBox::SET;
                    // Once tagged, no need to check other radii for this cell
                    break;
                }
            }
        }
    }
};

#endif /* EXTRACTIONTAGGER_HPP_ */
