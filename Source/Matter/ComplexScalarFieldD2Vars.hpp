/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef COMPLEXSCALARFIELDD2VARS_HPP_
#define COMPLEXSCALARFIELDD2VARS_HPP_

#include "CCZ4D2Vars.hpp"
#include "FourthOrderDerivatives.hpp"
#include "StateVariables.hpp"
#include "Tensor.hpp"

//! Second derivatives of the two real scalar components.
class ComplexScalarFieldD2Vars : public CCZ4D2Vars
{
  public:
    AMREX_GPU_DEVICE
    ComplexScalarFieldD2Vars(int ix, int iy, int iz,
                             const amrex::Array4<const amrex::Real> &state,
                             const FourthOrderDerivatives &a_deriv)
        : CCZ4D2Vars(ix, iy, iz, state, a_deriv)
    {
        phi1 = a_deriv.diff2_scalar(ix, iy, iz, state, c_phi);
        phi2 = a_deriv.diff2_scalar(ix, iy, iz, state, c_phi2);
    }

    Tensor<2, amrex::Real> phi1;
    Tensor<2, amrex::Real> phi2;
};

#endif /* COMPLEXSCALARFIELDD2VARS_HPP_ */
