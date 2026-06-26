/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef CCZ4D2VARS_HPP_
#define CCZ4D2VARS_HPP_

#include "FourthOrderDerivatives.hpp"
#include "SixthOrderDerivatives.hpp"
#include "StateVariables.hpp"
#include "Tensor.hpp"
#include "AMReX_Array4.H"

template <class deriv_t>
class CCZ4D2Vars
{
  public:
    AMREX_GPU_DEVICE CCZ4D2Vars(int ix, int iy, int iz,
                                const amrex::Array4<const amrex::Real> &state,
                                const deriv_t &a_deriv)
    {
        // Calculate the d2 quantities for all vars
        calculate_d2_derivs(ix, iy, iz, state, a_deriv);
    }

    // default empty contructor - for tests
    AMREX_GPU_HOST_DEVICE CCZ4D2Vars()
    {
        FOR (k, l)
        {
            chi[k][l]   = 0.0;
            lapse[k][l] = 0.0;
            FOR (i)
            {
                shift[i][k][l] = 0.0;
                FOR (j)
                {
                    h[i][j][k][l] = 0.0;
                }
            }
        }
    }

    AMREX_GPU_DEVICE void
    calculate_d2_derivs(int ix, int iy, int iz,
                        const amrex::Array4<const amrex::Real> &state,
                        const deriv_t &a_deriv)
    {
        // Calculate the d2 quantities for all required vars to calculate rhs
        chi   = a_deriv.diff2_scalar(ix, iy, iz, state, c_chi);
        lapse = a_deriv.diff2_scalar(ix, iy, iz, state, c_lapse);
        shift = a_deriv.diff2_vector(ix, iy, iz, state, c_shift1);
        h     = a_deriv.diff2_tensor(ix, iy, iz, state, c_h11);
    }

    Tensor<4, amrex::Real> h;
    Tensor<3, amrex::Real> shift;
    Tensor<2, amrex::Real> chi;
    Tensor<2, amrex::Real> lapse;
};

#endif /* CCZ4D2VARS_HPP */
