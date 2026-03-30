/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef CHITAGGER_HPP_
#define CHITAGGER_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "SphericalExtraction.hpp"
#include "Tensor.hpp"

#include <AMReX_Array4.H>
#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>
#include <AMReX_TagBox.H>

class ChiTagger
{
  protected:
    double m_dx;
    FourthOrderDerivatives m_deriv;
    amrex::Real m_threshold;

  public:
    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
    ChiTagger(const double dx, const amrex::Real a_threshold)
        : m_dx(dx), m_deriv(dx), m_threshold(a_threshold)
    {
    }

    AMREX_GPU_DEVICE void
    operator()(int i, int j, int k,
               const amrex::Array4<amrex::TagBox::TagType> &tags,
               const amrex::Array4<amrex::Real const> &state) const
    {
        const auto d2_chi = m_deriv.diff2_array_scalar(i, j, k, state, c_chi);
        amrex::Real mod_d2_chi = 0;
        FOR (idir, jdir)
        {
            mod_d2_chi += d2_chi(idir, jdir) * d2_chi(idir, jdir);
        }
        amrex::Real criterion = m_dx * std::sqrt(mod_d2_chi);

        if (criterion >= m_threshold)
        {
            tags(i, j, k) = amrex::TagBox::SET;
        }
    }
};

#endif /* CHITAGGER_HPP_ */
