/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SCALARFIELDD2VARS_HPP_
#define SCALARFIELDD2VARS_HPP_

#include "FourthOrderDerivatives.hpp"
#include "ScalarFieldVars.hpp"
#include "StateVariables.hpp"
#include "Tensor.hpp"

template <class deriv_t> class ScalarFieldD2Vars
{
  public:
    AMREX_GPU_DEVICE
    ScalarFieldD2Vars(int ix, int iy, int iz,
                      const amrex::Array4<const amrex::Real> &state,
                      const deriv_t &a_deriv)
    {
        // Calculate the d2 quantities for all vars needed for RHS
        phi = a_deriv.diff2_scalar(ix, iy, iz, state, c_phi);
    }

    TensorArray::Rank1Sym phi{};
};

#endif /* SCALARFIELDD2VARS_HPP */
