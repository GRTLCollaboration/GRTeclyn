/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(KLEINGORDONRHS_HPP_)
#error "This file should only be included through KleinGordonRHS.hpp"
#endif

#ifndef KLEINGORDONRHS_IMPL_HPP_
#define KLEINGORDONRHS_IMPL_HPP_

#include "KleinGordonRHS.hpp"

template <class model_t, class deriv_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
KleinGordonRHS<model_t, deriv_t>::operator()(
    int ix, int iy, int iz, const amrex::Array4<amrex::Real> &rhs_state,
    const amrex::Array4<amrex::Real const> &state) const
{
    const auto *state_ptr_ijk = state.ptr(ix, iy, iz);
    amrex::Array1D<amrex::Real, 0, AMREX_SPACEDIM>
        d2_phi{}; // no cross second order derivatives needed
    amrex::Array1D<int, 0, AMREX_SPACEDIM> strides{
        AMREX_D_DECL(1, static_cast<int>(state.stride.a[0]),
                     static_cast<int>(state.stride.a[1]))};

    FOR (i)
    {
        d2_phi(i) = m_deriv.diff2(state_ptr_ijk + c_phi * state.stride.a[2],
                                  strides(i));
    }

    rhs_equation(rhs_state.cellData(ix, iy, iz), state.cellData(ix, iy, iz),
                 d2_phi);

    m_deriv.add_dissipation(ix, iy, iz, rhs_state.cellData(ix, iy, iz), state,
                            m_sigma);
}

template <class model_t, class deriv_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
KleinGordonRHS<model_t, deriv_t>::rhs_equation(
    const amrex::CellData<amrex::Real> &rhs_cell_data,
    const amrex::CellData<amrex::Real const> &state_cell_data,
    const amrex::Array1D<amrex::Real, 0, AMREX_SPACEDIM> &d2_phi) const
{
    rhs_cell_data[c_phi] = state_cell_data[c_Pi];

    rhs_cell_data[c_Pi] = d2_phi.sum();

    // add on the potential
    amrex::Real V_of_phi = 0.0;
    amrex::Real dVdphi   = 0.0;

    m_model.compute_potential(V_of_phi, dVdphi, state_cell_data[c_phi]);

    rhs_cell_data[c_Pi] += dVdphi;
}

#endif // KLEINGORDONRHS_IMPL_HPP_