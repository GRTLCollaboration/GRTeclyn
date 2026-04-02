/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef CHIEXTRACTIONTAGGER_HPP_
#define CHIEXTRACTIONTAGGER_HPP_

#include "CCZ4Vars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "SphericalExtraction.hpp"

//! This class tags cells based on two criteria - the
//! value of the second derivs and the extraction regions
class ChiExtractionTagger
{
  protected:
    double m_dx;
    FourthOrderDerivatives m_deriv;
    amrex::Real m_threshold;
    // const SphericalExtraction::params_t m_params;  not GPU friendly
    int m_num_extraction_radii;
    const double *m_extraction_radii_ptr; // need to access via pointers on GPU
    const int *m_extraction_levels_ptr;   // need to access via pointers on GPU
    std::array<double, AMREX_SPACEDIM> m_center;
    int m_level;

  public:

    // The constructor
    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    ChiExtractionTagger(const double dx, const int a_level,
                        const amrex::Real a_threshold,
                        const SphericalExtraction::params_t &a_params,
                        const bool activate_extraction = false)
        : m_dx(dx), m_deriv(dx), m_threshold(a_threshold),
          m_num_extraction_radii(a_params.num_extraction_radii()),
          m_extraction_radii_ptr(a_params.extraction_radii().data()),
          m_extraction_levels_ptr(a_params.extraction_levels.data()),
          m_center(a_params.center), m_level(a_level)
    // NOLINTEND(bugprone-easily-swappable-parameters)
    {
        if (!activate_extraction)
        {
            // To avoid conditionals in the kernel
            m_num_extraction_radii = 0;
        }
    }

    AMREX_GPU_DEVICE void
    operator()(int ix, int iy, int iz,
               const amrex::Array4<amrex::TagBox::TagType> &tags,
               const amrex::Array4<amrex::Real const> &state) const
    {
        // first test the gradients for regions of high curvature
        const TensorArray::Rank1Sym d2_chi =
            m_deriv.diff2_sym_scalar(ix, iy, iz, state, c_chi);
        amrex::Real mod_d2_chi = 0;
        for (int i = 0; i < d2_chi.len(); i++)
        {
            mod_d2_chi += d2_chi(i) * d2_chi(i);
        }
        amrex::Real criterion = m_dx * std::sqrt(mod_d2_chi);
        if (criterion >= m_threshold)
        {
            tags(ix, iy, iz) = amrex::TagBox::SET;
        }

        // if extracting weyl data at a given radius, enforce a given resolution
        // there
        amrex::IntVect cell(AMREX_D_DECL(ix, iy, iz));
        for (int iradius = 0; iradius < m_num_extraction_radii; ++iradius)
        {
            // regrid if within extraction level and not at required
            // refinement
            if (m_level < m_extraction_levels_ptr[iradius])
            {
                const Coordinates coords(cell, m_dx, m_center);
                const amrex::Real r = coords.get_radius();
                // add a 20% buffer to extraction zone so not too near to
                // boundary
                if (r < 1.2 * m_extraction_radii_ptr[iradius])
                {
                    tags(ix, iy, iz) = amrex::TagBox::SET;
                }
            }
        }
    }
};

#endif /* CHIEXTRACTIONTAGGER_HPP_ */
