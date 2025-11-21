/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef CCZ4D2VARS_HPP_
#define CCZ4D2VARS_HPP_

#include "FourthOrderDerivatives.hpp"
#include "StateVariables.hpp"
#include "Tensor.hpp"
#include "AMReX_Array4.H"

class CCZ4D2Vars
{
  public:
    AMREX_GPU_DEVICE CCZ4D2Vars(int ix, int iy, int iz,
                                const amrex::Array4<const amrex::Real> &state,
                                const FourthOrderDerivatives &a_deriv)
    {
        // Calculate the d2 quantities for all vars
        calculate_d2_derivs(ix, iy, iz, state, a_deriv);
    }

    // default empty contructor - for tests
    AMREX_GPU_HOST_DEVICE CCZ4D2Vars()
    {
        FOR (k, l)
        {
            chi[k][l] = 0.0;
            // Theta[k][l] = 0.0;
            // K[k][l] = 0.0;
            lapse[k][l] = 0.0;
            FOR (i)
            {
                shift[i][k][l] = 0.0;
                // B[i][k][l] = 0.0;
                // Gamma[i][k][l] = 0.0;
                FOR (j)
                {
                    h[i][j][k][l] = 0.0;
                    // A[i][j][k][l] = 0.0;
                }
            }
        }
    }

    AMREX_GPU_DEVICE void
    calculate_d2_derivs(int ix, int iy, int iz,
                        const amrex::Array4<const amrex::Real> &state,
                        const FourthOrderDerivatives &a_deriv)
    {
        // Calculate the d2 quantities for all required vars to calculate rhs
        chi = a_deriv.diff2_scalar(ix, iy, iz, state, c_chi);
        // K     = a_deriv.diff2_scalar(ix, iy, iz, state, c_K);
        lapse = a_deriv.diff2_scalar(ix, iy, iz, state, c_lapse);
        // Theta = a_deriv.diff2_scalar(ix, iy, iz, state, c_Theta);
        shift = a_deriv.diff2_vector(ix, iy, iz, state, c_shift1);
        // Gamma = a_deriv.diff2_vector(ix, iy, iz, state, c_Gamma1);
        // B     = a_deriv.diff2_vector(ix, iy, iz, state, c_B1);
        h = a_deriv.diff2_tensor(ix, iy, iz, state, c_h11);
        // A     = a_deriv.diff2_tensor(ix, iy, iz, state, c_A11);
    }

    Tensor<4, amrex::Real> h;
    // Tensor<4, amrex::Real> A;
    // Tensor<3, amrex::Real> Gamma;
    Tensor<3, amrex::Real> shift;
    // Tensor<3, amrex::Real> B;
    Tensor<2, amrex::Real> chi;
    // Tensor<2, amrex::Real> K;
    Tensor<2, amrex::Real> lapse;
    // Tensor<2, amrex::Real> Theta;
};

#endif /* CCZ4D2VARS_HPP */