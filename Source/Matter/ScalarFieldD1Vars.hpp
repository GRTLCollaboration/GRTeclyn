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
    // NOLINTBEGIN(cppcoreguidelines-pro-type-member-init)
    AMREX_GPU_DEVICE
    ScalarFieldD1Vars(int ix, int iy, int iz,
                      const amrex::Array4<const amrex::Real> &state,
                      const FourthOrderDerivatives &a_deriv)
        : CCZ4D1Vars(ix, iy, iz, state, a_deriv)
    {
        // Calculate the d1 quantities for all vars
        auto d1_phi = a_deriv.diff1_array_scalar(ix, iy, iz, state, c_phi);
        auto d1_Pi  = a_deriv.diff1_array_scalar(ix, iy, iz, state, c_Pi);

        FOR (idir)
        {
            m_scalar_d1_state(c_phi - NUM_CCZ4_VARS, idir) = d1_phi(idir);
            m_scalar_d1_state(c_Pi - NUM_CCZ4_VARS, idir)  = d1_Pi(idir);
        }
    }
    // NOLINTEND(cppcoreguidelines-pro-type-member-init)

    // There are two scalar variables
    amrex::Array2D<amrex::Real, 0, 2, 0, 3> m_scalar_d1_state;

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &phi(int i) const
    {
        return m_scalar_d1_state(c_phi - NUM_CCZ4_VARS, i);
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &Pi(int i) const
    {
        return m_scalar_d1_state(c_Pi - NUM_CCZ4_VARS, i);
    }
};

#endif /* SCALARFIELDD1VARS_HPP */
