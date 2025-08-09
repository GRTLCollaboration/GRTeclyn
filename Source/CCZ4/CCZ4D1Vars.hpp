/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef CCZ4D1VARS_HPP_
#define CCZ4D1VARS_HPP_

#include "CCZ4Vars2.hpp"
#include "FourthOrderDerivatives2.hpp"
#include "StateVariables.hpp"
#include "Tensor.hpp"

class CCZ4D1Vars
{
  public:
    AMREX_GPU_DEVICE inline CCZ4D1Vars(
        int ix, int iy, int iz, const amrex::Array4<const amrex::Real> &state,
        const FourthOrderDerivatives2 &a_deriv)
    {
        // Calculate the d1 quantities for all vars
        chi   = a_deriv.diff1_scalar(ix, iy, iz, state, c_chi);
        K     = a_deriv.diff1_scalar(ix, iy, iz, state, c_K);
        lapse = a_deriv.diff1_scalar(ix, iy, iz, state, c_lapse);
        Theta = a_deriv.diff1_scalar(ix, iy, iz, state, c_Theta);
        shift = a_deriv.diff1_vector(ix, iy, iz, state, c_shift1);
        Gamma = a_deriv.diff1_vector(ix, iy, iz, state, c_Gamma1);
        B     = a_deriv.diff1_vector(ix, iy, iz, state, c_B1);
        h     = a_deriv.diff1_tensor(ix, iy, iz, state, c_h11);
        A     = a_deriv.diff1_tensor(ix, iy, iz, state, c_A11);
    }

    Tensor<3, amrex::Real> h;
    Tensor<3, amrex::Real> A;
    Tensor<2, amrex::Real> Gamma;
    Tensor<2, amrex::Real> shift;
    Tensor<2, amrex::Real> B;
    Tensor<1, amrex::Real> chi;
    Tensor<1, amrex::Real> K;
    Tensor<1, amrex::Real> lapse;
    Tensor<1, amrex::Real> Theta;
};

#endif /* CCZ4D1VARS_HPP */