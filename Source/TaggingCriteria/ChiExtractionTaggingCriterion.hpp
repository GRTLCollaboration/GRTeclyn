/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CHIEXTRACTIONTAGGINGCRITERION_HPP_
#define CHIEXTRACTIONTAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "SphericalExtraction.hpp"
#include "Tensor.hpp"

//! This class tags cells based on two criteria - the
//! value of the second derivs and the extraction regions
template <int NMAX = 2> // xxxxx TODO: max number of extractions
class ChiExtractionTaggingCriterion
{
  protected:
    const double m_dx;
    const FourthOrderDerivatives m_deriv;
    // const SphericalExtraction::params_t m_params;  not GPU friendly
    const int m_num_extraction_radii;
    std::array<double, NMAX> m_extraction_radii;
    std::array<int, NMAX> m_extraction_levels;
    const std::array<double, AMREX_SPACEDIM> m_center;
    const int m_level;
    const bool m_activate_extraction;

  public:
    template <class data_t> struct Vars
    {
        data_t chi; //!< Conformal factor

        template <typename mapping_function_t>
        AMREX_GPU_DEVICE void enum_mapping(mapping_function_t mapping_function)
        {
            using namespace VarsTools; // define_enum_mapping is part of
                                       // VarsTools
            define_enum_mapping(mapping_function, c_chi, chi);
        }
    };

    // The constructor
    ChiExtractionTaggingCriterion(const double dx, const int a_level,
                                  const SphericalExtraction::params_t a_params,
                                  const bool activate_extraction = false)
        : m_dx(dx), m_deriv(dx),
          m_num_extraction_radii(a_params.num_extraction_radii),
          m_center(a_params.center), m_level(a_level),
          m_activate_extraction(activate_extraction)
    {
        AMREX_ALWAYS_ASSERT(m_num_extraction_radii <= NMAX);
        for (int i = 0; i < m_num_extraction_radii; ++i)
        {
            m_extraction_radii[i]  = a_params.extraction_radii[i];
            m_extraction_levels[i] = a_params.extraction_levels[i];
        }
    }

    template <class data_t>
    AMREX_GPU_DEVICE data_t operator()(
        int i, int j, int k, amrex::Array4<data_t const> const &state) const
    {
        // first test the gradients for regions of high curvature
        const auto d2     = m_deriv.template diff2<Vars>(i, j, k, state);
        data_t mod_d2_chi = 0;
        FOR (idir, jdir)
        {
            mod_d2_chi += d2.chi[idir][jdir] * d2.chi[idir][jdir];
        }
        data_t criterion = m_dx * std::sqrt(mod_d2_chi);

        // if extracting weyl data at a given radius, enforce a given resolution
        // there
        if (m_activate_extraction)
        {
            amrex::IntVect cell(AMREX_D_DECL(i, j, k));
            for (int iradius = 0; iradius < m_num_extraction_radii; ++iradius)
            {
                // regrid if within extraction level and not at required
                // refinement
                if (m_level < m_extraction_levels[iradius])
                {
                    const Coordinates<data_t> coords(cell, m_dx, m_center);
                    const data_t r = coords.get_radius();
                    // add a 20% buffer to extraction zone so not too near to
                    // boundary
                    if (r < 1.2 * m_extraction_radii[iradius])
                    {
                        criterion = 100.0;
                    }
                }
            }
        }

        return criterion;
    }
};

#endif /* CHIEXTRACTIONTAGGINGCRITERION_HPP_ */
