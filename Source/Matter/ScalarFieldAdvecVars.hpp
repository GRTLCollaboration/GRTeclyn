/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SCALARFIELDADVECVARS_HPP_
#define SCALARFIELDADVECVARS_HPP_

#include "CCZ4AdvecVars.hpp"
#include "FourthOrderDerivatives.hpp"
#include "StateVariables.hpp"
#include "Tensor.hpp"

class ScalarFieldAdvecVars : public CCZ4AdvecVars
{
  public:
    AMREX_GPU_DEVICE inline ScalarFieldAdvecVars(
        int ix, int iy, int iz, const amrex::Array4<const amrex::Real> &state,
        const FourthOrderDerivatives &a_deriv)
        : CCZ4AdvecVars(ix, iy, iz, state, a_deriv)
    {
        Tensor<1, amrex::Real> shift_vector;
        FOR (idir)
        {
            shift_vector[idir] = state(ix, iy, iz, c_shift1 + idir);
        }

        // Calculate the d1 quantities for all vars
        phi = a_deriv.advec_scalar(ix, iy, iz, state, shift_vector, c_phi);
        Pi  = a_deriv.advec_scalar(ix, iy, iz, state, shift_vector, c_Pi);
    }

    amrex::Real phi;
    amrex::Real Pi;
};

#endif /* SCALARFIELDADVECVARS_HPP */