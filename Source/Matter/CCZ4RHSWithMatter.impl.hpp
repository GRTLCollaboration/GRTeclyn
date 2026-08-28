/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(CCZ4RHSWITHMATTER_HPP_)
#error "This file should only be included through MatterCCZ4RHS.hpp"
#endif

#ifndef CCZ4RHSWITHMATTER_IMPL_HPP_
#define CCZ4RHSWITHMATTER_IMPL_HPP_
#include "DimensionDefinitions.hpp"

template <class matter_t, class deriv_t>
CCZ4RHSWithMatter<matter_t, deriv_t>::CCZ4RHSWithMatter(amrex::Real a_dx)
    : CCZ4RHS<deriv_t>(a_dx, 0.0 /*No cosmological constant*/)
{
}

// Function to add in EM Tensor matter terms to CCZ4 rhs
template <class matter_t, class deriv_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
CCZ4RHSWithMatter<matter_t, deriv_t>::add_emtensor_rhs(
    const int ix, const int iy, const int iz,
    const amrex::Array4<amrex::Real> &rhs_state,
    const amrex::Array4<const amrex::Real> &state) const
{
    const amrex::CellData<amrex::Real> &rhs_cell_data =
        rhs_state.cellData(ix, iy, iz);
    const amrex::CellData<const amrex::Real> &state_cell_data =
        state.cellData(ix, iy, iz);

    const typename matter_t::Vars vars(state_cell_data);

    const auto h_UU = CCZ4Geometry::compute_inverse_metric(vars);

    const auto source = m_matter.compute_einstein_sources(ix, iy, iz, state,
                                                          this->m_deriv, h_UU);

    // Select the matter source terms without branching in the GPU kernel.
    const amrex::Real ccz4_coeff = 1.0 - this->m_params.bssn_coeff;

    const amrex::Real ccz4_K_matter_rhs =
        0.5 * vars.lapse() * (source.trS - 3.0 * source.rho);
    const amrex::Real bssn_K_matter_rhs =
        0.5 * vars.lapse() * (source.trS + source.rho);
    rhs_cell_data[c_K] += ccz4_coeff * ccz4_K_matter_rhs +
                          this->m_params.bssn_coeff * bssn_K_matter_rhs;

    const amrex::Real ccz4_Theta_matter_rhs = -vars.lapse() * source.rho;
    rhs_cell_data[c_Theta] =
        ccz4_coeff * (rhs_cell_data[c_Theta] + ccz4_Theta_matter_rhs);

    // Update RHS for other variables
    Tensor::Rank2 S_TF = source.S;

    CCZ4Geometry::make_trace_free(S_TF, vars, h_UU);

    FOR2_SYM(i, j)
    {

        rhs_cell_data[sym_var_idx(c_A11, i, j)] -=
            vars.chi() * vars.lapse() * S_TF(i, j);
    }

    FOR (i)
    {
        amrex::Real matter_term_Gamma = 0.0;
        FOR (j)
        {
            matter_term_Gamma -= 2.0 * vars.lapse() * h_UU(i, j) * source.j(j);
        }
        rhs_cell_data[c_Gamma1 + i] += matter_term_Gamma;
    }
}

template <class matter_t, class deriv_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
CCZ4RHSWithMatter<matter_t, deriv_t>::add_matter_rhs(
    int ix, int iy, int iz, const amrex::Array4<amrex::Real> &rhs_state,
    const amrex::Array4<const amrex::Real> &state) const
{
    m_matter.add_matter_rhs(ix, iy, iz, rhs_state, state, this->m_deriv);
}

template <class matter_t, class deriv_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
CCZ4RHSWithMatter<matter_t, deriv_t>::apply_dissipation(
    int ix, int iy, int iz, const amrex::Array4<amrex::Real> &rhs_state,
    const amrex::Array4<const amrex::Real> &state) const
{
    const amrex::CellData<amrex::Real> &rhs_cell_data =
        rhs_state.cellData(ix, iy, iz);
    this->m_deriv.add_dissipation(ix, iy, iz, rhs_cell_data, state,
                                  this->m_sigma, NUM_VARS);
}

#endif /* CCZ4RHSWITHMATTER_IMPL_HPP_ */
