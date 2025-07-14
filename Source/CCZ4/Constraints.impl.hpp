/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(CONSTRAINTS_HPP_)
#error "This file should only be included through Constraints.hpp"
#endif

#ifndef CONSTRAINTS_IMPL_HPP_
#define CONSTRAINTS_IMPL_HPP_

#include "Constraints.hpp"
#include "DimensionDefinitions.hpp"
#include "GRInterval.hpp"
#include "VarsTools.hpp"

// AMReX includes
#include <AMReX_AmrLevel.H>

template <class vars_t, class d1_vars_t, class d2_vars_t>
AMREX_GPU_DEVICE Constraints::Vars Constraints::constraint_equations(
    const vars_t &vars, const d1_vars_t &d1, const d2_vars_t &d2,
    const Tensor<2, amrex::Real> &h_UU, const chris_t &chris) const
{
    Vars out;

    if (m_c_Ham >= 0 || m_c_Ham_abs_terms >= 0)
    {
        auto ricci = CCZ4Geometry::compute_ricci(vars, d1, d2, h_UU, chris);

        auto A_UU         = TensorAlgebra::raise_all(vars.A, h_UU);
        amrex::Real tr_A2 = TensorAlgebra::compute_trace(vars.A, A_UU);

        out.Ham = ricci.scalar +
                  (GR_SPACEDIM - 1.) * vars.K * vars.K / GR_SPACEDIM - tr_A2;
        out.Ham -= 2 * m_cosmological_constant;

        out.Ham_abs_terms =
            std::abs(ricci.scalar) + std::abs(tr_A2) +
            std::abs((GR_SPACEDIM - 1.) * vars.K * vars.K / GR_SPACEDIM);
        out.Ham_abs_terms += 2.0 * std::abs(m_cosmological_constant);
    }

    if (m_c_Moms.size() > 0 || m_c_Moms_abs_terms.size() > 0)
    {
        Tensor<3, amrex::Real> covd_A;
        FOR (i, j, k)
        {
            covd_A[i][j][k] = d1.A[j][k][i];
            FOR (l)
            {
                covd_A[i][j][k] += -chris.ULL[l][i][j] * vars.A[l][k] -
                                   chris.ULL[l][i][k] * vars.A[l][j];
            }
        }
        FOR (i)
        {
            out.Mom[i]           = -(GR_SPACEDIM - 1.) * d1.K[i] / GR_SPACEDIM;
            out.Mom_abs_terms[i] = std::abs(out.Mom[i]);
        }
        Tensor<1, amrex::Real> covd_A_term = 0.0;
        Tensor<1, amrex::Real> d1_chi_term = 0.0;
        FOR (i, j, k)
        {
            covd_A_term[i] += h_UU[j][k] * covd_A[k][j][i];
            d1_chi_term[i] += -GR_SPACEDIM * h_UU[j][k] * vars.A[i][j] *
                              d1.chi[k] / (2 * vars.chi);
        }
        FOR (i)
        {
            out.Mom[i] += covd_A_term[i] + d1_chi_term[i];
            out.Mom_abs_terms[i] +=
                std::abs(covd_A_term[i]) + std::abs(d1_chi_term[i]);
        }
    }
    return out;
}

#endif /* CONSTRAINTS_IMPL_HPP_ */
