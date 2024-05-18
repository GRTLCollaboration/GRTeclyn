/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(MATTERCONSTRAINTS_HPP_)
#error "This file should only be included through MatterConstraints.hpp"
#endif

#ifndef MATTERCONSTRAINTS_IMPL_HPP_
#define MATTERCONSTRAINTS_IMPL_HPP_
#include "DimensionDefinitions.hpp"

template <class matter_t>
MatterConstraints<matter_t>::MatterConstraints(const matter_t a_matter,
                                               double dx, double G_Newton)
    : Constraints(dx, 0.0 /*No cosmological constant*/), my_matter(a_matter),
      m_G_Newton(G_Newton)
{
}

template <class matter_t>
template <class data_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
MatterConstraints<matter_t>::compute(int i, int j, int k,
                                     amrex::Array4<data_t> &cst,
                                     amrex::Array4<data_t const> &state) const
{
    // Load local vars and calculate derivs
    const auto vars = load_vars<BSSNMatterVars>(state.CellData(i, j, k));
    const auto d1 =
        m_deriv.template diff1<BSSNMatterVars>(state.CellData(i, j, k));
    const auto d2 =
        m_deriv.template diff2<BSSNMatterVars>(state.CellData(i, j, k));

    // Get the non matter terms for the constraints
    Vars<data_t> out = constraint_equations(vars, d1, d2);

    // Inverse metric and Christoffel symbol
    const auto h_UU  = TensorAlgebra::compute_inverse_sym(vars.h);
    const auto chris = TensorAlgebra::compute_christoffel(d1.h, h_UU);

    // Energy Momentum Tensor
    const auto emtensor = my_matter.compute_emtensor(vars, d1, h_UU, chris.ULL);

    // Hamiltonian constraint
    out.Ham += -16.0 * M_PI * m_G_Newton * emtensor.rho;

    // Momentum constraints
    FOR (i)
    {
        out.Mom[i] += -8.0 * M_PI * m_G_Newton * emtensor.Si[i];
    }

    // Write the constraints into the output FArrayBox
    store_vars(out, cst.cellData(i, j, k));
}

#endif /* MATTERCONSTRAINTS_IMPL_HPP_ */
