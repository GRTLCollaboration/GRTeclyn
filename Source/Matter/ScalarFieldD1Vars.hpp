/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SCALARFIELDD1VARS_HPP_
#define SCALARFIELDD1VARS_HPP_

#include "CCZ4D1Vars.hpp"
#include "FourthOrderDerivatives.hpp"
#include "ScalarFieldVars.hpp"
#include "StateVariables.hpp"
#include "Tensor.hpp"

class ScalarFieldD1Vars : public CCZ4D1Vars
{
  public:
    AMREX_GPU_DEVICE
    ScalarFieldD1Vars(int ix, int iy, int iz,
                      const amrex::Array4<const amrex::Real> &state,
                      const FourthOrderDerivatives &a_deriv)
        : CCZ4D1Vars(ix, iy, iz, state, a_deriv)
    {
        // Calculate the d1 quantities for all vars
        m_d1_state = a_deriv.diff1_state(ix, iy, iz, state);
    }

    amrex::GpuArray<Tensor<1, amrex::Real>, NUM_VARS> m_d1_state;

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &phi(int i) const
    {
        return m_d1_state[c_phi][i];
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &Pi(int i) const
    {
        return m_d1_state[c_Pi][i];
    }
};

#endif /* SCALARFIELDD1VARS_HPP */