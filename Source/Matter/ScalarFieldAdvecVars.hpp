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
    // NOLINTBEGIN(cppcoreguidelines-pro-type-member-init)
    AMREX_GPU_DEVICE
    ScalarFieldAdvecVars(int ix, int iy, int iz,
                         const amrex::Array4<const amrex::Real> &state,
                         const FourthOrderDerivatives &a_deriv)
        : CCZ4AdvecVars(ix, iy, iz, state, a_deriv)
    {
        // Calculate the advec quantities for all vars
        m_scalar_advec_state(c_phi - NUM_CCZ4_VARS) =
            a_deriv.advec_scalar(ix, iy, iz, state, m_shift_vector, c_phi);
        m_scalar_advec_state(c_Pi - NUM_CCZ4_VARS) =
            a_deriv.advec_scalar(ix, iy, iz, state, m_shift_vector, c_Pi);
    }
    // NOLINTEND(cppcoreguidelines-pro-type-member-init)

    // There are two scalar variables
    amrex::Array1D<amrex::Real, 0, 2> m_scalar_advec_state;

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &phi() const
    {
        return m_scalar_advec_state(c_phi - NUM_CCZ4_VARS);
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &Pi() const
    {
        return m_scalar_advec_state(c_Pi - NUM_CCZ4_VARS);
    }
};

#endif /* SCALARFIELDADVECVARS_HPP */
