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

template <class matter_t, class gauge_t, class deriv_t>
CCZ4RHSWithMatter<matter_t, gauge_t, deriv_t>::CCZ4RHSWithMatter(
    CCZ4_params_t<typename gauge_t::params_t> a_params, double a_dx,
    double a_sigma, int a_formulation, double a_G_Newton)
    : CCZ4RHS<gauge_t, deriv_t>(a_params, a_dx, a_sigma, a_formulation,
                                0.0 /*No cosmological constant*/),
      m_matter(matter_t()), m_G_Newton(a_G_Newton)
{
}

template <class matter_t, class gauge_t, class deriv_t>
CCZ4RHSWithMatter<matter_t, gauge_t, deriv_t>::CCZ4RHSWithMatter(
    matter_t a_matter, CCZ4_params_t<typename gauge_t::params_t> a_params, double a_dx,
    double a_sigma, int a_formulation, double a_G_Newton)
    : CCZ4RHS<gauge_t, deriv_t>(a_params, a_dx, a_sigma, a_formulation,
                                0.0 /*No cosmological constant*/),
      m_matter(a_matter), m_G_Newton(a_G_Newton)
{
}

template <class matter_t, class gauge_t, class deriv_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
CCZ4RHSWithMatter<matter_t, gauge_t, deriv_t>::operator()(
    int ix, int iy, int iz, const amrex::Array4<amrex::Real> &rhs_state,
    const amrex::Array4<amrex::Real const> &state) const
{

    const amrex::CellData<const amrex::Real> &state_cell_data =
        state.cellData(ix, iy, iz);
    const typename matter_t::Vars vars(state_cell_data);

    // Get the derivatives
    const typename matter_t::D1Vars d1(ix, iy, iz, state, this->m_deriv);
    const typename matter_t::D2Vars d2(ix, iy, iz, state, this->m_deriv);
    const typename matter_t::AdvecVars advec(ix, iy, iz, state, this->m_deriv);

    const amrex::CellData<amrex::Real> &rhs_cell_data =
        rhs_state.cellData(ix, iy, iz);

    this->rhs_equation(rhs_cell_data, vars, d1, d2, advec);

    // add RHS matter terms from EM Tensor
    add_emtensor_rhs(rhs_cell_data, vars, d1);

    // add evolution of matter fields themselves
    m_matter.add_matter_rhs(rhs_cell_data, vars, d1, d2, advec);

    // Add dissipation to all terms
    this->m_deriv.add_dissipation(ix, iy, iz, rhs_cell_data, state,
                                  this->m_sigma);
}

// Function to add in EM Tensor matter terms to CCZ4 rhs
template <class matter_t, class gauge_t, class deriv_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
CCZ4RHSWithMatter<matter_t, gauge_t, deriv_t>::add_emtensor_rhs(
    const amrex::CellData<amrex::Real> &rhs,
    const typename matter_t::Vars &vars,
    const typename matter_t::D1Vars &d1) const
{
    const auto h_UU  = CCZ4Geometry::compute_inverse_metric(vars);
    const auto chris = CCZ4Geometry::compute_christoffel(d1, h_UU);

    // Calculate elements of the decomposed stress energy tensor
    const auto emtensor = m_matter.compute_emtensor(vars, d1, h_UU, chris.ULL);

    // Update RHS for K and Theta depending on formulation
    if (this->m_formulation == CCZ4RHS<>::USE_BSSN)
    {
        rhs[c_K] += 4.0 * M_PI * m_G_Newton * vars.lapse() *
                    (emtensor.trS + emtensor.rho);
        rhs[c_Theta] = 0.0;
    }
    else
    {
        rhs[c_K] += 4.0 * M_PI * m_G_Newton * vars.lapse() *
                    (emtensor.trS - 3 * emtensor.rho);
        rhs[c_Theta] += -8.0 * M_PI * m_G_Newton * vars.lapse() * emtensor.rho;
    }

    // Update RHS for other variables
    Tensor<2, amrex::Real> S_TF = emtensor.S;
    CCZ4Geometry::make_trace_free(S_TF, vars, h_UU);

    FOR2_SYM(i, j)
    {
        rhs[VAR_IDX(c_A11, i, j)] +=
            -8.0 * M_PI * m_G_Newton * vars.chi() * vars.lapse() * S_TF[i][j];
    }

    FOR (i)
    {
        amrex::Real matter_term_Gamma = 0.0;
        FOR (j)
        {
            matter_term_Gamma += -16.0 * M_PI * m_G_Newton * vars.lapse() *
                                 h_UU[i][j] * emtensor.j[j];
        }
        rhs[c_Gamma1 + i] += matter_term_Gamma;
    }

    // Add matter contribution to RHS of gauge evolution
    this->m_gauge.rhs_gauge_add_matter_terms(rhs, vars, h_UU, emtensor,
                                             m_G_Newton);
}

#endif /* CCZ4RHSWITHMATTER_IMPL_HPP_ */
