/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SCALARFIELDADVECVARS_HPP_
#define SCALARFIELDADVECVARS_HPP_

#include "FourthOrderDerivatives.hpp"
#include "StateVariables.hpp"
#include "Tensor.hpp"

template <class deriv_t> class ScalarFieldAdvecVars
{
  public:
    // NOLINTBEGIN(cppcoreguidelines-pro-type-member-init)
    AMREX_GPU_DEVICE
    ScalarFieldAdvecVars(int ix, int iy, int iz,
                         const amrex::Array4<const amrex::Real> &state,
                         const TensorArray::Rank1 shift_vector,
                         const deriv_t &a_deriv)
    {
        // Calculate the advec quantities for all vars
        phi = a_deriv.advec_scalar(ix, iy, iz, state, shift_vector, c_phi);
        Pi  = a_deriv.advec_scalar(ix, iy, iz, state, shift_vector, c_Pi);
    }
    // NOLINTEND(cppcoreguidelines-pro-type-member-init)

    // There are two scalar variables
    amrex::Real phi{0.0};
    amrex::Real Pi{0.0};
};

#endif /* SCALARFIELDADVECVARS_HPP */
