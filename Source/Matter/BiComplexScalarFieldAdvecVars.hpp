/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef BICOMPLEXSCALARFIELDADVECVARS_HPP_
#define BICOMPLEXSCALARFIELDADVECVARS_HPP_

#include "CCZ4AdvecVars.hpp"
#include "FourthOrderDerivatives.hpp"
#include "StateVariables.hpp"
#include "Tensor.hpp"

//! Advection (shift-dotted) derivatives of the eight real scalar components.
class BiComplexScalarFieldAdvecVars : public CCZ4AdvecVars
{
  public:
    // NOLINTBEGIN(cppcoreguidelines-pro-type-member-init)
    AMREX_GPU_DEVICE
    BiComplexScalarFieldAdvecVars(
        int ix, int iy, int iz,
        const amrex::Array4<const amrex::Real> &state,
        const FourthOrderDerivatives &a_deriv)
        : CCZ4AdvecVars(ix, iy, iz, state, a_deriv)
    {
        m_phi1p = a_deriv.advec_scalar(ix, iy, iz, state, m_shift_vector, c_phi);
        m_Pi1p  = a_deriv.advec_scalar(ix, iy, iz, state, m_shift_vector, c_Pi);
        m_phi2p =
            a_deriv.advec_scalar(ix, iy, iz, state, m_shift_vector, c_phi2);
        m_Pi2p  = a_deriv.advec_scalar(ix, iy, iz, state, m_shift_vector, c_Pi2);
        m_phi1m =
            a_deriv.advec_scalar(ix, iy, iz, state, m_shift_vector, c_phi_m);
        m_Pi1m =
            a_deriv.advec_scalar(ix, iy, iz, state, m_shift_vector, c_Pi_m);
        m_phi2m =
            a_deriv.advec_scalar(ix, iy, iz, state, m_shift_vector, c_phi2_m);
        m_Pi2m =
            a_deriv.advec_scalar(ix, iy, iz, state, m_shift_vector, c_Pi2_m);
    }
    // NOLINTEND(cppcoreguidelines-pro-type-member-init)

    amrex::Real m_phi1p, m_Pi1p, m_phi2p, m_Pi2p;
    amrex::Real m_phi1m, m_Pi1m, m_phi2m, m_Pi2m;

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &
    phi1p() const { return m_phi1p; }
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &
    Pi1p() const { return m_Pi1p; }
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &
    phi2p() const { return m_phi2p; }
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &
    Pi2p() const { return m_Pi2p; }
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &
    phi1m() const { return m_phi1m; }
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &
    Pi1m() const { return m_Pi1m; }
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &
    phi2m() const { return m_phi2m; }
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &
    Pi2m() const { return m_Pi2m; }
};

#endif /* BICOMPLEXSCALARFIELDADVECVARS_HPP_ */
