/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef COMPLEXSCALARFIELDD1VARS_HPP_
#define COMPLEXSCALARFIELDD1VARS_HPP_

#include "CCZ4D1Vars.hpp"
#include "FourthOrderDerivatives.hpp"
#include "StateVariables.hpp"
#include "Tensor.hpp"

//! First derivatives of the two real scalar components.
class ComplexScalarFieldD1Vars : public CCZ4D1Vars
{
  public:
    // NOLINTBEGIN(cppcoreguidelines-pro-type-member-init)
    AMREX_GPU_DEVICE
    ComplexScalarFieldD1Vars(int ix, int iy, int iz,
                             const amrex::Array4<const amrex::Real> &state,
                             const FourthOrderDerivatives &a_deriv)
        : CCZ4D1Vars(ix, iy, iz, state, a_deriv)
    {
        m_phi1_d1 = a_deriv.diff1_scalar(ix, iy, iz, state, c_phi);
        m_phi2_d1 = a_deriv.diff1_scalar(ix, iy, iz, state, c_phi2);
    }
    // NOLINTEND(cppcoreguidelines-pro-type-member-init)

    Tensor<1, amrex::Real> m_phi1_d1;
    Tensor<1, amrex::Real> m_phi2_d1;

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &phi1(int i) const
    {
        return m_phi1_d1[i];
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &phi2(int i) const
    {
        return m_phi2_d1[i];
    }
};

#endif /* COMPLEXSCALARFIELDD1VARS_HPP_ */
