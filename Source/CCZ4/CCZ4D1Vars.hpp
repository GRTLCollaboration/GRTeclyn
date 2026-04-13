/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef CCZ4D1VARS_HPP_
#define CCZ4D1VARS_HPP_

#include "FourthOrderDerivatives.hpp"
#include "StateVariables.hpp"
#include "Tensor.hpp"
#include "AMReX_Array4.H"

class CCZ4D1Vars
{
  public:
    // NOLINTBEGIN(cppcoreguidelines-pro-type-member-init)
    AMREX_GPU_DEVICE CCZ4D1Vars(int ix, int iy, int iz,
                                const amrex::Array4<const amrex::Real> &state,
                                const FourthOrderDerivatives &a_deriv)
    {
        m_d1_state = a_deriv.diff1_state<NUM_CCZ4_VARS>(ix, iy, iz, state);
    }

    // empty constructor used for tests
    AMREX_GPU_HOST_DEVICE CCZ4D1Vars()
    {
        for (int ivar = 0; ivar < NUM_CCZ4_VARS; ++ivar)
        {
            FOR (i)
            {
                m_d1_state(ivar, i) = 0.0;
            }
        }
    }
    // NOLINTEND(cppcoreguidelines-pro-type-member-init)

    amrex::Array2D<amrex::Real, 0, NUM_CCZ4_VARS, 0, AMREX_SPACEDIM>
        m_d1_state{};

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &chi(int i) const
    {
        return m_d1_state(c_chi, i);
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &h(int i, int j,
                                                             int k) const
    {
        return m_d1_state(VAR_IDX(c_h11, i, j), k);
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &K(int i) const
    {
        return m_d1_state(c_K, i);
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &A(int i, int j,
                                                             int k) const
    {
        return m_d1_state(VAR_IDX(c_A11, i, j), k);
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &Theta(int i) const
    {
        return m_d1_state(c_Theta, i);
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &Gamma(int i,
                                                                 int j) const
    {
        return m_d1_state(c_Gamma1 + i, j);
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &lapse(int i) const
    {
        return m_d1_state(c_lapse, i);
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &shift(int i,
                                                                 int j) const
    {
        return m_d1_state(c_shift1 + i, j);
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &B(int i, int j) const
    {
        return m_d1_state(c_B1 + i, j);
    }
};

#endif /* CCZ4D1VARS_HPP */
