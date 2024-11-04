/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(KLEINGORDONRHS_HPP_)
#error "This file should only be included through KleinGordonRHS.hpp"
#endif

#ifndef KLEINGORDONRHS_IMPL_HPP_
#define KLEINGORDONRHS_IMPL_HPP_

#include "KleinGordonRHS.hpp"

template <class deriv_t>
template <class data_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
KleinGordonRHS<deriv_t>::compute(int i, int j, int k,
                                 const amrex::Array4<data_t const> &input,
                                 const amrex::Array4<data_t> &output) const

{
    const auto vars = load_vars<Vars>(input.cellData(i, j, k));
    const auto d1   = this->m_deriv.template diff1<Vars>(i, j, k, input);
    const auto d2   = this->m_deriv.template diff2<Diff2Vars>(i, j, k, input);

    Vars<data_t> rhs;
    rhs_equation(rhs, vars, d1, d2);

    m_deriv.add_dissipation(i, j, k, rhs, input, m_sigma);

    store_vars(output.cellData(i, j, k), rhs);
}

template <class deriv_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void KleinGordonRHS<deriv_t>::rhs_equation(
    vars_t<data_t> &rhs, const vars_t<data_t> &vars,
    const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2) const
{
    rhs.phi = vars.Pi;
    rhs.Pi  = TensorAlgebra::compute_trace(d2.phi) - std::sin(vars.phi);
}

#endif
