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

        for (int i = 0; i < TensorArray::Rank1Sym::len(); i++)
        {
            chi(i)   = 0.;
            lapse(i) = 0.;

            FOR (k)
            {
                shift(k, i) = 0.;
            }


            for (int j = 0; j < TensorArray::Rank2Sym::ylen(); j++)
            {
                h(i, j) = 0.;
            }
        }
    }

    AMREX_GPU_DEVICE void
    calculate_d2_derivs(int ix, int iy, int iz,
                        const amrex::Array4<const amrex::Real> &state,
                        const FourthOrderDerivatives &a_deriv)
    {
        // Calculate the d2 quantities for all required vars to calculate rhs
        chi   = a_deriv.diff2_scalar(ix, iy, iz, state, c_chi);
        lapse = a_deriv.diff2_scalar(ix, iy, iz, state, c_lapse);
        shift = a_deriv.diff2_vector(ix, iy, iz, state, c_shift1);
        h     = a_deriv.diff2_tensor(ix, iy, iz, state, c_h11);
    }

    TensorArray::Rank2Sym h;
    amrex::Array2D<amrex::Real, 0, 3, 0, 6> shift;
    TensorArray::Rank1Sym chi;
    TensorArray::Rank1Sym lapse;
};

#endif /* CCZ4D2VARS_HPP */
