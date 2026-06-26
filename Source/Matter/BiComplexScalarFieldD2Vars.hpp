/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef BICOMPLEXSCALARFIELDD2VARS_HPP_
#define BICOMPLEXSCALARFIELDD2VARS_HPP_

#include "CCZ4D2Vars.hpp"
#include "FourthOrderDerivatives.hpp"
#include "StateVariables.hpp"
#include "Tensor.hpp"

//! Second derivatives of the four real scalar components.
class BiComplexScalarFieldD2Vars : public CCZ4D2Vars
{
  public:
    AMREX_GPU_DEVICE
    BiComplexScalarFieldD2Vars(int ix, int iy, int iz,
                               const amrex::Array4<const amrex::Real> &state,
                               const FourthOrderDerivatives &a_deriv)
        : CCZ4D2Vars(ix, iy, iz, state, a_deriv)
    {
        phi1p = a_deriv.diff2_scalar(ix, iy, iz, state, c_phi);
        phi2p = a_deriv.diff2_scalar(ix, iy, iz, state, c_phi2);
        phi1m = a_deriv.diff2_scalar(ix, iy, iz, state, c_phi_m);
        phi2m = a_deriv.diff2_scalar(ix, iy, iz, state, c_phi2_m);
    }

    Tensor<2, amrex::Real> phi1p;
    Tensor<2, amrex::Real> phi2p;
    Tensor<2, amrex::Real> phi1m;
    Tensor<2, amrex::Real> phi2m;
};

#endif /* BICOMPLEXSCALARFIELDD2VARS_HPP_ */
