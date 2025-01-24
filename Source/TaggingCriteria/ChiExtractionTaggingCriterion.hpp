/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
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
class ChiExtractionTaggingCriterion
{
  protected:
    double m_dx;
    FourthOrderDerivatives m_deriv;
    // const SphericalExtraction::params_t m_params;  not GPU friendly
    int m_num_extraction_radii;
    const double *m_extraction_radii_ptr; // need to access via pointers on GPU
    const int *m_extraction_levels_ptr;   // need to access via pointers on GPU
    std::array<double, AMREX_SPACEDIM> m_center;
    int m_level;

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
    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    ChiExtractionTaggingCriterion(const double dx, const int a_level,
                                  const SphericalExtraction::params_t &a_params,
                                  const bool activate_extraction = false)
        : m_dx(dx), m_deriv(dx), m_center(a_params.center), m_level(a_level)
    // NOLINTEND(bugprone-easily-swappable-parameters)
    {
        // To avoid conditionals in the kernel, we set the number of extraction
        // radii to 0 if activate_extraction == false
        m_num_extraction_radii =
            (activate_extraction) ? a_params.num_extraction_radii() : 0;
        m_extraction_radii_ptr  = a_params.extraction_radii().data();
        m_extraction_levels_ptr = a_params.extraction_levels.data();
    }

    template <class data_t>
    AMREX_GPU_DEVICE data_t operator()(
        int i, int j, int k, const amrex::Array4<data_t const> &state) const
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
        amrex::IntVect cell(AMREX_D_DECL(i, j, k));
        for (int iradius = 0; iradius < m_num_extraction_radii; ++iradius)
        {
            // regrid if within extraction level and not at required
            // refinement
            if (m_level < m_extraction_levels_ptr[iradius])
            {
                const Coordinates<data_t> coords(cell, m_dx, m_center);
                const data_t r = coords.get_radius();
                // add a 20% buffer to extraction zone so not too near to
                // boundary
                if (r < 1.2 * m_extraction_radii_ptr[iradius])
                {
                    criterion = 100.0;
                }
            }
        }

        return criterion;
    }
};

#endif /* CHIEXTRACTIONTAGGINGCRITERION_HPP_ */
