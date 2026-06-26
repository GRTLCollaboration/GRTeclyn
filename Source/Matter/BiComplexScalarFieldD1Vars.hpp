/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef BICOMPLEXSCALARFIELDD1VARS_HPP_
#define BICOMPLEXSCALARFIELDD1VARS_HPP_

#include "CCZ4D1Vars.hpp"
#include "FourthOrderDerivatives.hpp"
#include "StateVariables.hpp"
#include "Tensor.hpp"

//! First derivatives of the four real scalar components (phi1p, phi2p of the
//! canonical field and phi1m, phi2m of the phantom field).
class BiComplexScalarFieldD1Vars : public CCZ4D1Vars
{
  public:
    // NOLINTBEGIN(cppcoreguidelines-pro-type-member-init)
    AMREX_GPU_DEVICE
    BiComplexScalarFieldD1Vars(int ix, int iy, int iz,
                               const amrex::Array4<const amrex::Real> &state,
                               const FourthOrderDerivatives &a_deriv)
        : CCZ4D1Vars(ix, iy, iz, state, a_deriv)
    {
        m_phi1p_d1 = a_deriv.diff1_scalar(ix, iy, iz, state, c_phi);
        m_phi2p_d1 = a_deriv.diff1_scalar(ix, iy, iz, state, c_phi2);
        m_phi1m_d1 = a_deriv.diff1_scalar(ix, iy, iz, state, c_phi_m);
        m_phi2m_d1 = a_deriv.diff1_scalar(ix, iy, iz, state, c_phi2_m);
    }
    // NOLINTEND(cppcoreguidelines-pro-type-member-init)

    Tensor<1, amrex::Real> m_phi1p_d1;
    Tensor<1, amrex::Real> m_phi2p_d1;
    Tensor<1, amrex::Real> m_phi1m_d1;
    Tensor<1, amrex::Real> m_phi2m_d1;

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &
    phi1p(int i) const { return m_phi1p_d1[i]; }
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &
    phi2p(int i) const { return m_phi2p_d1[i]; }
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &
    phi1m(int i) const { return m_phi1m_d1[i]; }
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &
    phi2m(int i) const { return m_phi2m_d1[i]; }
};

#endif /* BICOMPLEXSCALARFIELDD1VARS_HPP_ */
