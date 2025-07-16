/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#include "Constraints.hpp"
#include "DimensionDefinitions.hpp"
#include "GRInterval.hpp"
#include "VarsTools.hpp"

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
Constraints::Constraints(double dx, int a_c_Ham, const Interval &a_c_Moms,
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
        "must calculate one of either Ham or Ham_abs_terms");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        a_c_Moms.size() * a_c_Moms_abs_terms.size() <= 0,
        "must choose at most one of Mom or Mom_abs_terms");
    const auto &moms_interval =
        (a_c_Moms.size() > 0) ? a_c_Moms : a_c_Moms_abs_terms;
    if (moms_interval.size() > 0)
    {
        AMREX_ALWAYS_ASSERT(moms_interval.size() == (s_calc_mom_norm ? 1 : 3));
    }
}
// NOLINTEND(bugprone-easily-swappable-parameters)

AMREX_GPU_DEVICE void
Constraints::compute(int i, int j, int k, const amrex::Array4<amrex::Real> &cst,
                     const amrex::Array4<amrex::Real const> &state) const
{
    const auto d1 = m_deriv.template diff1<MetricVars>(i, j, k, state);
    const auto d2 = m_deriv.template diff2<Diff2Vars>(i, j, k, state);

    const auto state_cell = state.cellData(i, j, k);
    const auto vars       = load_vars<MetricVars>(state_cell);
    const auto h_UU       = TensorAlgebra::compute_inverse_sym(vars.h);
    const auto chris      = TensorAlgebra::compute_christoffel(d1.h, h_UU);

    Vars out = constraint_equations(vars, d1, d2, h_UU, chris);

    const auto cst_cell = cst.cellData(i, j, k);
    store_vars(out, cst_cell);
}

AMREX_GPU_DEVICE void
Constraints::store_vars(const Vars &out,
                        const amrex::CellData<amrex::Real> &current_cell) const
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
        amrex::Real Mom_sq = 0.0;
        FOR (i)
        {
            Mom_sq += out.Mom[i] * out.Mom[i];
        }
        amrex::Real Mom                = sqrt(Mom_sq);
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
        amrex::Real Mom_abs_terms_sq = 0.0;
        FOR (i)
        {
            Mom_abs_terms_sq += out.Mom_abs_terms[i] * out.Mom_abs_terms[i];
        }
        amrex::Real Mom_abs_terms                = sqrt(Mom_abs_terms_sq);
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
        static_cast<int>(comp_names.size()), comp_names, compute_mf,
        [=](const amrex::Box &box) { return amrex::grow(box, num_ghosts); },
        &amrex::cell_quartic_interp);

    // We only need the non-gauge CCZ4 variables to calculate the constraints
    derive_lst.addComponent(name, desc_lst, a_state_index, 0, c_lapse);
}

void Constraints::compute_mf(amrex::MultiFab &out_mf, int dcomp, int ncomp,
                             const amrex::MultiFab &src_mf,
                             const amrex::Geometry &geomdata,
                             amrex::Real /*time*/, const int * /*bcrec*/,
                             int /*level*/)
{
    const auto &out_arrays = out_mf.arrays();
    const auto &src_arrays = src_mf.const_arrays();
    int iham               = dcomp;
    Interval imom          = Interval(dcomp + 1, dcomp + AMREX_SPACEDIM);
    Constraints constraints(geomdata.CellSize(0), iham, imom);

    amrex::ParallelFor(
        out_mf, out_mf.nGrowVect(),
        [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
        {
            constraints.compute(i, j, k, out_arrays[box_no],
                                src_arrays[box_no]);
        });
}
