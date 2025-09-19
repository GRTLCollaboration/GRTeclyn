/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(CONSTRAINTSWITHMATTER_HPP_)
#error "This file should only be included through ConstraintsWithMatter.hpp"
#endif

#ifndef CONSTRAINTSWITHMATTER_IMPL_HPP_
#define CONSTRAINTSWITHMATTER_IMPL_HPP_

#include "ConstraintsWithMatter.hpp"
#include "DimensionDefinitions.hpp"
#include "GRParmParse.hpp"

template <class matter_t>
ConstraintsWithMatter<matter_t>::ConstraintsWithMatter(
    double dx, double G_Newton, int a_c_Ham, const Interval &a_c_Moms,
    int a_c_Ham_abs_terms /* defaulted*/,
    const Interval &a_c_Moms_abs_terms /*defaulted*/)
    : Constraints(dx, a_c_Ham, a_c_Moms, a_c_Ham_abs_terms, a_c_Moms_abs_terms,
                  0.0 /*No cosmological constant*/),
      m_G_Newton(G_Newton)
{
}

template <class matter_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
ConstraintsWithMatter<matter_t>::operator()(
    int ix, int iy, int iz, const amrex::Array4<amrex::Real> &constraints,
    const amrex::Array4<amrex::Real const> &state) const
{
    const amrex::CellData<const amrex::Real> &state_cell_data =
        state.cellData(ix, iy, iz);
    typename matter_t::Vars vars(state_cell_data);
    const typename matter_t::D1Vars d1(ix, iy, iz, state, m_deriv);
    const Tensor<2, amrex::Real> d2_chi =
        m_deriv.diff2(ix, iy, iz, state, c_chi);
    const Tensor<4, amrex::Real> d2_h =
        m_deriv.diff2_tensor(ix, iy, iz, state, c_h11);

    // Inverse metric and Christoffel symbol
    const auto h_UU  = CCZ4Geometry::compute_inverse_metric(vars);
    const auto chris = CCZ4Geometry::compute_christoffel(d1, h_UU);

    // Get the non matter terms for the constraints
    constraints_t out =
        constraint_equations(vars, d1, d2_chi, d2_h, h_UU, chris);

    // Energy Momentum Tensor
    const auto emtensor = my_matter.compute_emtensor(vars, d1, h_UU, chris.ULL);

    // Hamiltonian constraint
    if (m_c_Ham >= 0 || m_c_Ham_abs_terms >= 0)
    {
        out.Ham           += -16.0 * M_PI * m_G_Newton * emtensor.rho;
        out.Ham_abs_terms += 16.0 * M_PI * m_G_Newton * std::abs(emtensor.rho);
    }

    // Momentum constraints
    if (m_c_Moms.size() > 0 || m_c_Moms_abs_terms.size() > 0)
    {
        FOR (i)
        {
            out.Mom(i) += -8.0 * M_PI * m_G_Newton * emtensor.j(i);
            out.Mom_abs_terms(i) +=
                8.0 * M_PI * m_G_Newton * std::abs(emtensor.j(i));
        }
    }
    // Write the constraints into the output FArrayBox
    store_vars(out, constraints.cellData(ix, iy, iz));
}

template <class matter_t>
void ConstraintsWithMatter<matter_t>::set_up(int a_state_index,
                                             bool a_calc_mom_norm)
{

    int num_ghosts = 2;

    auto &derive_lst     = amrex::AmrLevel::get_derive_lst();
    const auto &desc_lst = amrex::AmrLevel::get_desc_lst();

    const auto &comp_names = (a_calc_mom_norm) ? Constraints::var_names_norm
                                               : Constraints::var_names;
    // Add Constraints to the derive list
    derive_lst.add(
        Constraints::name, amrex::IndexType::TheCellType(),
        static_cast<int>(comp_names.size()), comp_names,
        ConstraintsWithMatter::compute_mf, [=](const amrex::Box &box)
        { return amrex::grow(box, num_ghosts); }, &amrex::cell_quartic_interp);

    derive_lst.addComponent(Constraints::name, desc_lst, a_state_index, 0,
                            NUM_VARS);
}
template <class matter_t>
void ConstraintsWithMatter<matter_t>::compute_mf(
    amrex::MultiFab &out_mf, int dcomp, int ncomp,
    const amrex::MultiFab &src_mf, const amrex::Geometry &geomdata,
    amrex::Real /*time*/, const int * /*bcrec*/, int /*level*/)
{
    const auto &out_arrays = out_mf.arrays();
    const auto &src_arrays = src_mf.const_arrays();

    GRParmParse pp;
    amrex::Real G_Newton = 0;

    pp.get("G_Newton", G_Newton, 0);

    amrex::Real dx = geomdata.CellSize(0);
    int iham       = dcomp; // Ham
    Interval imom =
        Interval(dcomp + 1, dcomp + AMREX_SPACEDIM); // Mom1, Mom2, Mom3

    AMREX_ALWAYS_ASSERT(ncomp == (1 + AMREX_SPACEDIM));

    ConstraintsWithMatter<matter_t> my_matter_constraints(dx, G_Newton, iham,
                                                          imom);

    amrex::ParallelFor(
        out_mf,
        [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz) noexcept
        {
            my_matter_constraints(ix, iy, iz, out_arrays[box_no],
                                  src_arrays[box_no]);
        });
}

#endif /* CONSTRAINTSWITHMATTER_IMPL_HPP_ */
