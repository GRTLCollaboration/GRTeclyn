/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(CONSTRAINTS_HPP_)
#error "This file should only be included through Constraints.hpp"
#endif

#ifndef CONSTRAINTS_IMPL_HPP_
#define CONSTRAINTS_IMPL_HPP_

#include "DimensionDefinitions.hpp"
#include "GRInterval.hpp"
#include "VarsTools.hpp"

// AMReX includes
#include <AMReX_AmrLevel.H>

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
inline Constraints::Constraints(
    double dx, int a_c_Ham, const Interval &a_c_Moms,
    int a_c_Ham_abs_terms /*defaulted*/,
    const Interval &a_c_Moms_abs_terms /*defaulted*/,
    double cosmological_constant /*defaulted*/)
    : m_deriv(dx), m_c_Ham(a_c_Ham), m_c_Moms(a_c_Moms),
      m_c_Ham_abs_terms(a_c_Ham_abs_terms),
      m_c_Moms_abs_terms(a_c_Moms_abs_terms),
      m_cosmological_constant(cosmological_constant)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        (a_c_Ham >= 0 && a_c_Ham_abs_terms < 0) ||
            (a_c_Ham < 0 && a_c_Ham_abs_terms >= 0),
        must calculate one of either Ham or Ham_abs_terms);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        a_c_Moms.size() * a_c_Moms_abs_terms.size() <= 0,
        must choose at most one of Mom or Mom_abs_terms);
    const auto &moms_interval =
        (a_c_Moms.size() > 0) ? a_c_Moms : a_c_Moms_abs_terms;
    if (moms_interval.size() > 0)
    {
        AMREX_ALWAYS_ASSERT(moms_interval.size() == (s_calc_mom_norm ? 1 : 3));
    }
}
// NOLINTEND(bugprone-easily-swappable-parameters)

template <class data_t>
AMREX_GPU_DEVICE void
Constraints::compute(int i, int j, int k, const amrex::Array4<data_t> &cst,
                     const amrex::Array4<data_t const> &state) const
{
    const auto d1 = m_deriv.template diff1<MetricVars>(i, j, k, state);
    const auto d2 = m_deriv.template diff2<Diff2Vars>(i, j, k, state);

    const auto state_cell = state.cellData(i, j, k);
    const auto vars       = load_vars<MetricVars>(state_cell);
    const auto h_UU       = TensorAlgebra::compute_inverse_sym(vars.h);
    const auto chris      = TensorAlgebra::compute_christoffel(d1.h, h_UU);

    Vars<data_t> out = constraint_equations(vars, d1, d2, h_UU, chris);

    const auto cst_cell = cst.cellData(i, j, k);
    store_vars(out, cst_cell);
}

template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t>
AMREX_GPU_DEVICE Constraints::Vars<data_t> Constraints::constraint_equations(
    const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2, const Tensor<2, data_t> &h_UU,
    const chris_t<data_t> &chris) const
{
    Vars<data_t> out;

    if (m_c_Ham >= 0 || m_c_Ham_abs_terms >= 0)
    {
        auto ricci = CCZ4Geometry::compute_ricci(vars, d1, d2, h_UU, chris);

        auto A_UU    = TensorAlgebra::raise_all(vars.A, h_UU);
        data_t tr_A2 = TensorAlgebra::compute_trace(vars.A, A_UU);

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
        Tensor<3, data_t> covd_A;
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
        Tensor<1, data_t> covd_A_term = 0.0;
        Tensor<1, data_t> d1_chi_term = 0.0;
        const data_t chi_regularised  = simd_max(1e-6, vars.chi);
        FOR (i, j, k)
        {
            covd_A_term[i] += h_UU[j][k] * covd_A[k][j][i];
            d1_chi_term[i] += -GR_SPACEDIM * h_UU[j][k] * vars.A[i][j] *
                              d1.chi[k] / (2 * chi_regularised);
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

template <class data_t>
AMREX_GPU_DEVICE void
Constraints::store_vars(const Vars<data_t> &out,
                        const amrex::CellData<data_t> &current_cell) const
{
    if (m_c_Ham >= 0)
    {
        current_cell[m_c_Ham] = out.Ham;
    }
    if (m_c_Ham_abs_terms >= 0)
    {
        current_cell[m_c_Ham_abs_terms] = out.Ham_abs_terms;
    }
    if (m_c_Moms.size() == GR_SPACEDIM)
    {
        FOR (i)
        {
            int ivar           = m_c_Moms.begin() + i;
            current_cell[ivar] = out.Mom[i];
        }
    }
    else if (m_c_Moms.size() == 1)
    {
        data_t Mom_sq = 0.0;
        FOR (i)
        {
            Mom_sq += out.Mom[i] * out.Mom[i];
        }
        data_t Mom                     = sqrt(Mom_sq);
        current_cell[m_c_Moms.begin()] = Mom;
    }
    if (m_c_Moms_abs_terms.size() == GR_SPACEDIM)
    {
        FOR (i)
        {
            int ivar           = m_c_Moms_abs_terms.begin() + i;
            current_cell[ivar] = out.Mom_abs_terms[i];
        }
    }
    else if (m_c_Moms_abs_terms.size() == 1)
    {
        data_t Mom_abs_terms_sq = 0.0;
        FOR (i)
        {
            Mom_abs_terms_sq += out.Mom_abs_terms[i] * out.Mom_abs_terms[i];
        }
        data_t Mom_abs_terms                     = sqrt(Mom_abs_terms_sq);
        current_cell[m_c_Moms_abs_terms.begin()] = Mom_abs_terms;
    }
}

void Constraints::set_up(int a_state_index, bool a_calc_mom_norm)
{
    s_calc_mom_norm = a_calc_mom_norm;
    int num_ghosts  = 2; // no advection terms so only need 2 ghost cells

    const auto &comp_names = (s_calc_mom_norm) ? var_names_norm : var_names;
    auto &derive_lst       = amrex::AmrLevel::get_derive_lst();
    const auto &desc_lst   = amrex::AmrLevel::get_desc_lst();

    derive_lst.add(
        name, amrex::IndexType::TheCellType(),
        static_cast<int>(comp_names.size()), comp_names,
        amrex::DeriveFuncFab(), // null function because we won't use
                                // it.
        [=](const amrex::Box &box) { return amrex::grow(box, num_ghosts); },
        &amrex::cell_quartic_interp);

    // We only need the non-gauge CCZ4 variables to calculate the constraints
    derive_lst.addComponent(name, desc_lst, a_state_index, 0, c_lapse);
}

#endif /* CONSTRAINTS_IMPL_HPP_ */
