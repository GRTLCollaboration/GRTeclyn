/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(WEYL4WITHMATTER_HPP_)
#error "This file should only be included through Weyl4WithMatter.hpp"
#endif

#ifndef WEYL4WITHMATTER_IMPL_HPP_
#define WEYL4WITHMATTER_IMPL_HPP_

template <class matter_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void Weyl4WithMatter<matter_t>::operator()(
    int ix, int iy, int iz, const amrex::Array4<amrex::Real> &weyl_scalars,
    const amrex::Array4<amrex::Real const> &state) const
{
    const amrex::CellData<const amrex::Real> &state_cell_data =
        state.cellData(ix, iy, iz);
    const typename matter_t::Vars vars(state_cell_data);

    //    const typename matter_t::D1Vars d1(ix, iy, iz, state, m_deriv);
    // we only need d2 of chi and h
    const Tensor::Sym12Rank2 d2_chi =
        m_deriv.d2_scalar(ix, iy, iz, state, c_chi);
    const Tensor::Sym12Sym34Rank4 d2_h =
        m_deriv.d2_tensor(ix, iy, iz, state, c_h11);

    auto d1_h        = m_deriv.d1_sym_tensor(ix, iy, iz, state, c_h11);
    const auto h_UU  = CCZ4Geometry::compute_inverse_metric(vars);
    const auto chris = CCZ4Geometry::compute_christoffel(d1_h, h_UU);

    // Get the coordinates
    const Coordinates coords(amrex::IntVect{AMREX_D_DECL(ix, iy, iz)}, m_dx,
                             m_center);

    // Compute the spatial volume element
    const auto epsilon3_LUU = compute_epsilon3_LUU(vars, h_UU);

    // Compute the E and B fields
    // This needs d1 chi, K, h, A
    auto d1_chi   = m_deriv.d1_scalar(ix, iy, iz, state, c_chi);
    auto d1_Gamma = m_deriv.d1_vector(ix, iy, iz, state, c_Gamma1);
    auto d1_K     = m_deriv.d1_scalar(ix, iy, iz, state, c_K);
    auto d1_A     = m_deriv.d1_sym_tensor(ix, iy, iz, state, c_A11);

    EBFields_t ebfields =
        compute_EB_fields(vars, d1_chi, d1_Gamma, d1_h, d1_K, d1_A, d2_chi,
                          d2_h, epsilon3_LUU, h_UU, chris);

    // Add in matter terms to E and B fields

    add_matter_EB(ebfields, ix, iy, iz, state, epsilon3_LUU, h_UU, chris);

    // work out the Newman Penrose scalar
    weyl_scalar_t out = compute_Weyl4(ebfields, vars, h_UU, coords);

    // Write the rhs into the output FArrayBox
    weyl_scalars(ix, iy, iz, m_out_comp)     = out.Real;
    weyl_scalars(ix, iy, iz, m_out_comp + 1) = out.Im;
}

template <class matter_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
Weyl4WithMatter<matter_t>::add_matter_EB(
    EBFields_t &ebfields, const int ix, const int iy, const int iz,
    const amrex::Array4<const amrex::Real> &state,
    const Tensor::Rank3 &epsilon3_LUU, const Tensor::Rank2 &h_UU,
    const chris_t &chris) const
{
    const amrex::CellData<const amrex::Real> &state_cell_data =
        state.cellData(ix, iy, iz);

    const typename matter_t::Vars vars(state_cell_data);

    // Calculate decomposed energy momentum tensor components
    const auto emtensor =
        m_matter.compute_emtensor(ix, iy, iz, state, m_deriv, h_UU);

    Tensor::Rank2 S_TF = emtensor.S;
    CCZ4Geometry::make_trace_free(S_TF, vars, h_UU);

    // as we made the vacuum expression of Bij explictly symmetric and Eij
    // explictly trace-free, only Eij has matter terms
    FOR (i, j)
    {
        ebfields.E(i, j) += -4.0 * M_PI * m_G_Newton * S_TF(i, j);
    }
}

template <class matter_t>
void Weyl4WithMatter<matter_t>::set_up(int a_state_index)
{

    int num_ghosts = 2;

    auto &derive_lst     = amrex::AmrLevel::get_derive_lst();
    const auto &desc_lst = amrex::AmrLevel::get_desc_lst();

    // Add MatterWeyl4 to the derive list
    derive_lst.add(
        Weyl4::name, amrex::IndexType::TheCellType(),
        static_cast<int>(Weyl4::var_names.size()), Weyl4::var_names,
        Weyl4WithMatter::compute_mf, [=](const amrex::Box &box)
        { return amrex::grow(box, num_ghosts); }, &amrex::cell_quartic_interp);

    derive_lst.addComponent(Weyl4::name, desc_lst, a_state_index, 0, NUM_VARS);
}

template <class matter_t>
void Weyl4WithMatter<matter_t>::compute_mf(amrex::MultiFab &out_mf,
                                           int out_comp, int ncomp,
                                           const amrex::MultiFab &src_mf,
                                           const amrex::Geometry &geomdata,
                                           amrex::Real /*time*/,
                                           const int * /*bcrec*/, int /*level*/)
{
    const auto &out_arrays = out_mf.arrays();
    const auto &src_arrays = src_mf.const_arrays();

    GRParmParse pp;
    amrex::Real G_Newton = 0;

    pp.queryAdd("G_newton", G_Newton);

    Weyl4WithMatter<matter_t> my_weyl4_with_matter(
        geomdata.CellSize(0), out_comp, G_Newton);

    amrex::ParallelFor(
        out_mf,
        [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz) noexcept
        {
            my_weyl4_with_matter(ix, iy, iz, out_arrays[box_no],
                                 src_arrays[box_no]);
        });
}

#endif /* WEYL4WITHMATTER_IMPL_HPP_ */
