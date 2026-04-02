/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SCALARFIELDD2VARS_HPP_
#define SCALARFIELDD2VARS_HPP_

#include "CCZ4D2Vars.hpp"
#include "FourthOrderDerivatives.hpp"
#include "ScalarFieldVars.hpp"
#include "StateVariables.hpp"
#include "Tensor.hpp"

class ScalarFieldD2Vars : public CCZ4D2Vars
{
  public:
    AMREX_GPU_DEVICE
    ScalarFieldD2Vars(int ix, int iy, int iz,
                      const amrex::Array4<const amrex::Real> &state,
                      const FourthOrderDerivatives &a_deriv)
        : CCZ4D2Vars(ix, iy, iz, state, a_deriv)
    {
        // Calculate the d2 quantities for all vars needed for RHS
        phi = a_deriv.diff2_sym_scalar(ix, iy, iz, state, c_phi);
    }

    TensorArray::Rank1Sym phi{};
};

#endif /* SCALARFIELDD2VARS_HPP */
