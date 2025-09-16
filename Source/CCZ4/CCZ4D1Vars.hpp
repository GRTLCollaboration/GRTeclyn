/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef CCZ4D1VARS_HPP_
#define CCZ4D1VARS_HPP_

#include "FourthOrderDerivatives.hpp"
#include "StateVariables.hpp"
#include "Tensor.hpp"
#include "AMReX_Array4.H"

class CCZ4D1Vars
{
  public:
    AMREX_GPU_DEVICE CCZ4D1Vars(int ix, int iy, int iz,
                                const amrex::Array4<const amrex::Real> &state,
                                const FourthOrderDerivatives &a_deriv)
    {
        // Calculate the d1 quantities for all vars
        calculate_d1_derivs(ix, iy, iz, state, a_deriv);
    }

    // default empty contructor
    AMREX_GPU_DEVICE CCZ4D1Vars() { zero_d1_derivs(); }

    // default empty contructor
    AMREX_GPU_DEVICE void zero_d1_derivs()
    {
        FOR (k)
        {
            chi[k]   = 0.0;
            Theta[k] = 0.0;
            K[k]     = 0.0;
            lapse[k] = 0.0;
            FOR (i)
            {
                shift[i][k] = 0.0;
                B[i][k]     = 0.0;
                Gamma[i][k] = 0.0;
                FOR (j)
                {
                    h[i][j][k] = 0.0;
                    A[i][j][k] = 0.0;
                }
            }
        }
    }

    AMREX_GPU_DEVICE void
    calculate_d1_derivs(int ix, int iy, int iz,
                        const amrex::Array4<const amrex::Real> &state,
                        const FourthOrderDerivatives &a_deriv)
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