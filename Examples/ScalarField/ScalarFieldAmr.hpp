/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SCALARFIELDAMR_HPP_
#define SCALARFIELDAMR_HPP_

#include "GRAmr.hpp"
#include "ParticleInterpolator.hpp"

//! AMR hierarchy carrying the interpolators used for line extraction.
class ScalarFieldAmr : public GRAmr
{
  public:
    ParticleInterpolator<1> phi_interpolator;
    ParticleInterpolator<1> rho_interpolator;

    explicit ScalarFieldAmr(amrex::LevelBld *a_level_bld) : GRAmr(a_level_bld)
    {
    }
};

#endif /* SCALARFIELDAMR_HPP_ */
