/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef COMPLEXSCALARFIELDADVECVARS_HPP_
#define COMPLEXSCALARFIELDADVECVARS_HPP_

#include "CCZ4AdvecVars.hpp"
#include "FourthOrderDerivatives.hpp"
#include "StateVariables.hpp"
#include "Tensor.hpp"

//! Advection (shift-dotted) derivatives of the two real scalar components.
class ComplexScalarFieldAdvecVars : public CCZ4AdvecVars
{
  public:
    // NOLINTBEGIN(cppcoreguidelines-pro-type-member-init)
    AMREX_GPU_DEVICE
    ComplexScalarFieldAdvecVars(int ix, int iy, int iz,
                                const amrex::Array4<const amrex::Real> &state,
                                const FourthOrderDerivatives &a_deriv)
        : CCZ4AdvecVars(ix, iy, iz, state, a_deriv)
    {
        m_phi1_advec =
            a_deriv.advec_scalar(ix, iy, iz, state, m_shift_vector, c_phi);
        m_Pi1_advec =
            a_deriv.advec_scalar(ix, iy, iz, state, m_shift_vector, c_Pi);
        m_phi2_advec =
            a_deriv.advec_scalar(ix, iy, iz, state, m_shift_vector, c_phi2);
        m_Pi2_advec =
            a_deriv.advec_scalar(ix, iy, iz, state, m_shift_vector, c_Pi2);
    }
    // NOLINTEND(cppcoreguidelines-pro-type-member-init)

    amrex::Real m_phi1_advec;
    amrex::Real m_Pi1_advec;
    amrex::Real m_phi2_advec;
    amrex::Real m_Pi2_advec;

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &phi1() const
    {
        return m_phi1_advec;
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &Pi1() const
    {
        return m_Pi1_advec;
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &phi2() const
    {
        return m_phi2_advec;
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &Pi2() const
    {
        return m_Pi2_advec;
    }
};

#endif /* COMPLEXSCALARFIELDADVECVARS_HPP_ */
