/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef EXTRACTIONTAGGER_HPP_
#define EXTRACTIONTAGGER_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "SphericalExtraction.hpp"
#include "Tensor.hpp"

#include <AMReX_Array4.H>
#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>
#include <AMReX_TagBox.H>

class ExtractionTagger
{
  protected:
    double m_dx;
    int m_num_extraction_radii;
    const double *m_extraction_radii_ptr;
    const int *m_extraction_levels_ptr;
    std::array<double, AMREX_SPACEDIM> m_center;
    int m_level;

  public:
    // Constructor takes only what it needs for extraction tagging
    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
    ExtractionTagger(const double dx, const int a_level,
                     const spherical_extraction_params_t &a_params,
                     const bool activate_extraction = false)
        : m_dx(dx), m_num_extraction_radii(a_params.num_extraction_radii()),
          m_extraction_radii_ptr(a_params.extraction_radii().data()),
          m_extraction_levels_ptr(a_params.extraction_levels.data()),
          m_center(a_params.center), m_level(a_level)
    {
        if (!activate_extraction)
        {
            // Avoids conditionals in the kernel by setting num to 0
            m_num_extraction_radii = 0;
        }
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
                const int exponent  = std::min(m_extraction_levels_ptr[iradius] - m_level - 1, 0);
                const double factor = std::pow(1.5, exponent);

                // Add a 20% buffer to extraction zone so not too near boundary
                if (r < factor * 1.2 * m_extraction_radii_ptr[iradius])
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
