/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef CCZ4ADVECVARS_HPP_
#define CCZ4ADVECVARS_HPP_

#include "FourthOrderDerivatives.hpp"
#include "StateVariables.hpp"
#include "Tensor.hpp"
#include "AMReX_Array4.H"

class CCZ4AdvecVars
{
  public:
    // NOLINTBEGIN(cppcoreguidelines-pro-type-member-init)
    AMREX_GPU_DEVICE
    CCZ4AdvecVars(int ix, int iy, int iz,
                  const amrex::Array4<const amrex::Real> &state,
                  const FourthOrderDerivatives &a_deriv)
    {
        Tensor<1, amrex::Real> shift_vector;
        FOR (idir)
        {
            shift_vector[idir] = state(ix, iy, iz, c_shift1 + idir);
        }

        // Calculate the advec quantities for all vars
        m_advec_state = a_deriv.advec_state(ix, iy, iz, state, shift_vector);
    }
    // NOLINTEND(cppcoreguidelines-pro-type-member-init)

    amrex::GpuArray<amrex::Real, NUM_VARS> m_advec_state;

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &chi() const
    {
        return m_advec_state[c_chi];
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &h(int i, int j) const
    {
        return m_advec_state[var_idx(c_h11, i, j)];
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &K() const
    {
        return m_advec_state[c_K];
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &A(int i, int j) const
    {
        return m_advec_state[var_idx(c_A11, i, j)];
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &Theta() const
    {
        return m_advec_state[c_Theta];
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &Gamma(int i) const
    {
        return m_advec_state[c_Gamma1 + i];
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &lapse() const
    {
        return m_advec_state[c_lapse];
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &shift(int i) const
    {
        return m_advec_state[c_shift1 + i];
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &B(int i) const
    {
        return m_advec_state[c_B1 + i];
    }
};

#endif /* CCZ4ADVECVARS_HPP */