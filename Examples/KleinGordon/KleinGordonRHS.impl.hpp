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
template <class data_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
KleinGordonRHS<model_t, deriv_t>::compute(
    int i, int j, int k, const amrex::Array4<data_t const> &input,
    const amrex::Array4<data_t> &output) const

{
    const auto *input_ptr_ijk = input.ptr(i, j, k);
    amrex::Array1D<amrex::Real, 0, AMREX_SPACEDIM>
        d2phi{}; // no cross second order derivatives needed
    amrex::Array1D<int, 0, AMREX_SPACEDIM> strides{AMREX_D_DECL(
        1, static_cast<int>(input.jstride), static_cast<int>(input.kstride))};

    FOR (i)
    {
        d2phi(i) =
            m_deriv.diff2(input_ptr_ijk + c_phi * input.nstride, 0, strides(i));
    }

    rhs_equation(input.cellData(i, j, k), output.cellData(i, j, k), d2phi);

    // add dissipation term
    amrex::Real phi_dissipation = 0.0;
    amrex::Real Pi_dissipation  = 0.0;

    FOR (i)
    {

        phi_dissipation +=
            m_sigma * m_deriv.dissipation_term(
                          input_ptr_ijk + c_phi * input.nstride, 0, strides(i));

        Pi_dissipation +=
            m_sigma * m_deriv.dissipation_term(
                          input_ptr_ijk + c_Pi * input.nstride, 0, strides(i));
    }

    output.cellData(i, j, k)[c_phi] += phi_dissipation;
    output.cellData(i, j, k)[c_Pi]  += Pi_dissipation;
}

template <class model_t, class deriv_t>
template <class data_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
KleinGordonRHS<model_t, deriv_t>::rhs_equation(
    const amrex::CellData<data_t const> &input_cell_data,
    const amrex::CellData<data_t> &output_cell_data,
    const amrex::Array1D<data_t, 0, AMREX_SPACEDIM> &d2phi) const
{
    output_cell_data[c_phi] = input_cell_data[c_Pi];

    output_cell_data[c_Pi] = d2phi.sum();

    // add on the potential
    data_t V_of_phi = 0.0;
    data_t dVdphi   = 0.0;

    m_model.compute_potential(V_of_phi, dVdphi, input_cell_data[c_phi]);

    output_cell_data[c_Pi] += dVdphi;
}

#endif // KLEINGORDONRHS_IMPL_HPP_
