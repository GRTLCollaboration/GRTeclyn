/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SCALARFIELDD1VARS_HPP_
#define SCALARFIELDD1VARS_HPP_

#include "FourthOrderDerivatives.hpp"
#include "ScalarFieldVars.hpp"
#include "StateVariables.hpp"
#include "Tensor.hpp"

template <class deriv_t = FourthOrderDerivatives> class ScalarFieldD1Vars
{
  public:
    // NOLINTBEGIN(cppcoreguidelines-pro-type-member-init)
    AMREX_GPU_DEVICE
    ScalarFieldD1Vars(int ix, int iy, int iz,
                      const amrex::Array4<const amrex::Real> &state,
                      const deriv_t &a_deriv)
    {
        // Calculate the d1 quantities for all vars
        phi = a_deriv.diff1_scalar(ix, iy, iz, state, c_phi);
        Pi  = a_deriv.diff1_scalar(ix, iy, iz, state, c_Pi);
    }
    // NOLINTEND(cppcoreguidelines-pro-type-member-init)

    // There are two scalar variables
    TensorArray::Rank1 phi{0.0, 0.0, 0.0};
    TensorArray::Rank1 Pi{0.0, 0.0, 0.0};
};

#endif /* SCALARFIELDD1VARS_HPP */
