/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef CCZ4ADVECVARS_HPP_
#define CCZ4ADVECVARS_HPP_

#include "CCZ4Vars2.hpp"
#include "FourthOrderDerivatives2.hpp"
#include "StateVariables.hpp"
#include "Tensor.hpp"

class CCZ4AdvecVars
{
  public:
    AMREX_GPU_DEVICE inline CCZ4AdvecVars(
        int ix, int iy, int iz, const amrex::Array4<const amrex::Real> &state,
        const FourthOrderDerivatives2 &a_deriv)
    {
        Tensor<1, amrex::Real> shift_vector;
        FOR (idir)
        {
            shift_vector[idir] = state(ix, iy, iz, c_shift1 + idir);
        }

        // Calculate the d1 quantities for all vars
        chi   = a_deriv.advec_scalar(ix, iy, iz, state, shift_vector, c_chi);
        K     = a_deriv.advec_scalar(ix, iy, iz, state, shift_vector, c_K);
        lapse = a_deriv.advec_scalar(ix, iy, iz, state, shift_vector, c_lapse);
        Theta = a_deriv.advec_scalar(ix, iy, iz, state, shift_vector, c_Theta);
        shift = a_deriv.advec_vector(ix, iy, iz, state, shift_vector, c_shift1);
        Gamma = a_deriv.advec_vector(ix, iy, iz, state, shift_vector, c_Gamma1);
        B     = a_deriv.advec_vector(ix, iy, iz, state, shift_vector, c_B1);
        h     = a_deriv.advec_tensor(ix, iy, iz, state, shift_vector, c_h11);
        A     = a_deriv.advec_tensor(ix, iy, iz, state, shift_vector, c_A11);
    }

    Tensor<2, amrex::Real> h;
    Tensor<2, amrex::Real> A;
    Tensor<1, amrex::Real> Gamma;
    Tensor<1, amrex::Real> shift;
    Tensor<1, amrex::Real> B;
    amrex::Real chi;
    amrex::Real K;
    amrex::Real lapse;
    amrex::Real Theta;
};

#endif /* CCZ4ADVECVARS_HPP */