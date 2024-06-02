/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(MATTERWEYL4_HPP_)
#error "This file should only be included through MatterWeyl4.hpp"
#endif

#ifndef MATTERWEYL4_IMPL_HPP_
#define MATTERWEYL4_IMPL_HPP_

template <class matter_t>
template <class data_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
MatterWeyl4<matter_t>::compute(int i, int j, int k,
                               const amrex::Array4<data_t> &derive,
                               const amrex::Array4<data_t const> &state) const
{

    // copy data from chombo gridpoint into local variables
    const auto vars = load_vars<Vars>(state.cellData(i, j, k));
    const auto d1   = m_deriv.template diff1<Vars>(i, j, k, state);
    const auto d2   = m_deriv.template diff2<Diff2Vars>(i, j, k, state);

    // Get the coordinates
    amrex::IntVect cell_coords(AMREX_D_DECL(i, j, k));
    const Coordinates<data_t> coords(cell_coords, m_dx, m_center);

    // Compute the inverse metric
    using namespace TensorAlgebra;
    const auto h_UU  = compute_inverse_sym(vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);

    // Compute the spatial volume element
    const auto epsilon3_LUU = compute_epsilon3_LUU(vars, h_UU);

    // Compute the E and B fields
    EBFields_t<data_t> ebfields =
        compute_EB_fields(vars, d1, d2, epsilon3_LUU, h_UU, chris);

    // Add in matter terms to E and B fields
    add_matter_EB(ebfields, vars, d1, epsilon3_LUU, h_UU, chris);

    // work out the Newman Penrose scalar
    NPScalar_t<data_t> out =
        compute_Weyl4(ebfields, vars, d1, d2, h_UU, coords);

    // Write the rhs into the output FArrayBox
    derive(i, j, k, m_dcomp)     = out.Real;
    derive(i, j, k, m_dcomp + 1) = out.Im;
    //    current_cell.store_vars(out.Real, c_Weyl4_Re);
    //    current_cell.store_vars(out.Im, c_Weyl4_Im);
}

template <class matter_t>
template <class data_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void MatterWeyl4<matter_t>::add_matter_EB(
    EBFields_t<data_t> &ebfields, const Vars<data_t> &vars,
    const Vars<Tensor<1, data_t>> &d1, const Tensor<3, data_t> &epsilon3_LUU,
    const Tensor<2, data_t> &h_UU, const chris_t<data_t> &chris) const
{
    // Calculate decomposed energy momentum tensor components
    const auto emtensor = m_matter.compute_emtensor(vars, d1, h_UU, chris.ULL);

    Tensor<2, data_t> Sij_TF = emtensor.Sij;
    TensorAlgebra::make_trace_free(Sij_TF, vars.h, h_UU);

    // as we made the vacuum expression of Bij explictly symmetric and Eij
    // explictly trace-free, only Eij has matter terms
    FOR (i, j)
    {
        ebfields.E[i][j] += -4.0 * M_PI * m_G_Newton * Sij_TF[i][j];
    }
}

#endif /* MATTERWEYL4_IMPL_HPP_ */
