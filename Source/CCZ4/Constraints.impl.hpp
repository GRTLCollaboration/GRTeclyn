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

// AMReX includes
#include <AMReX_AmrLevel.H>

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
AMREX_FORCE_INLINE
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
Constraints::operator()(int ix, int iy, int iz,
                        const amrex::Array4<amrex::Real> &constraints,
                        const amrex::Array4<const amrex::Real> &state) const
{
    // We do not want to amend the cell data values, so use const CCZ4Vars
    const amrex::CellData<const amrex::Real> &state_cell_data =
        state.cellData(ix, iy, iz);
    CCZ4Vars vars(state_cell_data);

    // we need d1 chi, K, h, A... this just gets all of them
    auto d1_chi   = m_deriv.d1_scalar(ix, iy, iz, state, c_chi);
    auto d1_Gamma = m_deriv.d1_vector(ix, iy, iz, state, c_Gamma1);
    auto d1_K     = m_deriv.d1_scalar(ix, iy, iz, state, c_K);
    auto d1_A     = m_deriv.d1_sym_tensor(ix, iy, iz, state, c_A11);
    auto d1_h     = m_deriv.d1_sym_tensor(ix, iy, iz, state, c_h11);

    // we only need d2 of chi and h
    const Tensor::Sym12Rank2 d2_chi =
        m_deriv.d2_scalar(ix, iy, iz, state, c_chi);
    const Tensor::Sym12Sym34Rank4 d2_h =
        m_deriv.d2_tensor(ix, iy, iz, state, c_h11);

    const auto h_UU  = CCZ4Geometry::compute_inverse_metric(vars);
    const auto chris = CCZ4Geometry::compute_christoffel(d1_h, h_UU);

    constraints_t out = constraint_equations(vars, d1_chi, d1_Gamma, d1_h, d1_K,
                                             d1_A, d2_chi, d2_h, h_UU, chris);

    // TODO: Simplify this storing so less choice but more readable
    const auto constraint_cell_data = constraints.cellData(ix, iy, iz);
    store_vars(out, constraint_cell_data);
}

AMREX_GPU_DEVICE
Constraints::constraints_t Constraints::constraint_equations(
    const CCZ4Vars &vars, const Tensor::Rank1 &d1_chi,
    const Tensor::Rank2 &d1_Gamma, const Tensor::Sym12Rank3 &d1_h,
    const Tensor::Rank1 &d1_K, const Tensor::Sym12Rank3 &d1_A,
    const Tensor::Sym12Rank2 &d2_chi, const Tensor::Sym12Sym34Rank4 &d2_h,
    const Tensor::Rank2 &h_UU, const chris_t &chris) const
{
    constraints_t out;

    if (m_c_Ham >= 0 || m_c_Ham_abs_terms >= 0)
    {
        auto ricci = CCZ4Geometry::compute_ricci(vars, d1_chi, d1_Gamma, d1_h,
                                                 d2_chi, d2_h, h_UU, chris);

        // This is A_ij A^ij
        amrex::Real Aij_squared = CCZ4Geometry::compute_Aij_squared(vars, h_UU);

        out.Ham = ricci.scalar +
                  (GR_SPACEDIM - 1.) * vars.K() * vars.K() / GR_SPACEDIM -
                  Aij_squared - 2.0 * m_cosmological_constant;

        out.Ham_abs_terms =
            std::abs(ricci.scalar) +
            std::abs((GR_SPACEDIM - 1.) * vars.K() * vars.K() / GR_SPACEDIM) +
            std::abs(Aij_squared) + 2.0 * std::abs(m_cosmological_constant);
    }

    if (m_c_Moms.size() > 0 || m_c_Moms_abs_terms.size() > 0)
    {
        // Covariant derivative of \bar A_ij
        Tensor::Rank3 covd_A{};
        FOR (i, j, k)
        {
            covd_A(i, j, k) = d1_A(j, k, i);
            FOR (l)
            {
                covd_A(i, j, k) += -chris.ULL(l, i, j) * vars.A(l, k) -
                                   chris.ULL(l, i, k) * vars.A(l, j);
            }
        }
        FOR (i)
        {
            out.Mom(i)           = -(GR_SPACEDIM - 1.) * d1_K(i) / GR_SPACEDIM;
            out.Mom_abs_terms(i) = std::abs(out.Mom(i));
        }
        Tensor::Rank1 covd_A_term{};
        Tensor::Rank1 d1_chi_term{};
        FOR (i, j, k)
        {
            covd_A_term(i) += h_UU(j, k) * covd_A(k, j, i);
            d1_chi_term(i) += -GR_SPACEDIM * h_UU(j, k) * vars.A(i, j) *
                              d1_chi(k) / (2.0 * vars.chi());
        }
        FOR (i)
        {
            out.Mom(i) += covd_A_term(i) + d1_chi_term(i);
            out.Mom_abs_terms(i) +=
                std::abs(covd_A_term(i)) + std::abs(d1_chi_term(i));
        }
    }

    return out;
}

AMREX_GPU_DEVICE void
Constraints::store_vars(const constraints_t &out,
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
            current_cell[ivar] = out.Mom(i);
        }
    }
    else if (m_c_Moms.size() == 1)
    {
        amrex::Real Mom_sq = 0.0;
        FOR (i)
        {
            Mom_sq += out.Mom(i) * out.Mom(i);
        }
        amrex::Real Mom                = std::sqrt(Mom_sq);
        current_cell[m_c_Moms.begin()] = Mom;
    }
    if (m_c_Moms_abs_terms.size() == GR_SPACEDIM)
    {
        FOR (i)
        {
            int ivar           = m_c_Moms_abs_terms.begin() + i;
            current_cell[ivar] = out.Mom_abs_terms(i);
        }
    }
    else if (m_c_Moms_abs_terms.size() == 1)
    {
        amrex::Real Mom_abs_terms_sq = 0.0;
        FOR (i)
        {
            Mom_abs_terms_sq += out.Mom_abs_terms(i) * out.Mom_abs_terms(i);
        }
        amrex::Real Mom_abs_terms                = std::sqrt(Mom_abs_terms_sq);
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

    // Get all the vars to allow us to use the CCZ4Vars class
    // Technically not all needed but probably doesn't hurt performance
    derive_lst.addComponent(name, desc_lst, a_state_index, 0, NUM_CCZ4_VARS);
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
        [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz) noexcept
        { constraints(ix, iy, iz, out_arrays[box_no], src_arrays[box_no]); });
}

#endif /* CONSTRAINTS_IMPL_HPP_ */
