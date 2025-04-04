/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(CCZ4RHS_HPP_)
#error "This file should only be included through CCZ4RHS.hpp"
#endif

#ifndef CCZ4RHS_IMPL_HPP_
#define CCZ4RHS_IMPL_HPP_

#include "Macros.hpp"
#include "GRInterval.hpp"
#include "VarsTools.hpp"
#include "VarsWrapper.hpp"

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
template <class gauge_t, class deriv_t>
inline CCZ4RHS<gauge_t, deriv_t>::CCZ4RHS(
    CCZ4_params_t<typename gauge_t::params_t> a_params, double a_dx,
    double a_sigma, int a_formulation, double a_cosmological_constant)
    : m_params(a_params), m_gauge(a_params), m_sigma(a_sigma),
      m_formulation(a_formulation),
      m_cosmological_constant(a_cosmological_constant), m_deriv(a_dx)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    // A user who wants to use BSSN should also have damping paramters = 0
    if (m_formulation == USE_BSSN)
    {
        if ((m_params.kappa1 != 0.) || (m_params.kappa2 != 0.) ||
            (m_params.kappa3 != 0.))
        {
            amrex::Abort("BSSN formulation is selected - CCZ4 kappa values "
                         "should be set to zero in params");
        }
    }
    if (m_formulation > USE_BSSN)
    {
        amrex::Abort("The requested formulation is not supported");
    }
}

template <class gauge_t, class deriv_t>
template <class data_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CCZ4RHS<gauge_t, deriv_t>::compute(
    int i, int j, int k, const amrex::Array4<data_t> &rhs,
    const amrex::Array4<data_t const> &state) const
{
    const auto advec =
        m_deriv.template advection<data_t, NUM_CCZ4_VARS>(i, j, k, state);

    const auto d1   = m_deriv.template diff1<data_t, NUM_CCZ4_VARS>(i, j, k, state);

    Tensor<2, data_t> diff2_lapse = m_deriv.template diff2<data_t>(i, j, k, state, c_lapse);
    Tensor<2, data_t> diff2_chi = m_deriv.template diff2<data_t>(i, j, k, state, c_chi);

    Tensor<1, Tensor<2, data_t>> diff2_shift;
    FOR(idx)
    {
        diff2_shift[idx] = m_deriv.template diff2<data_t>(i, j, k, state, c_shift1 + idx);
    }

    Tensor<2, Tensor<2, data_t>> diff2_h;
    FOR(idx1, idx2)
    {
        diff2_h[idx1][idx2] = m_deriv.template diff2<data_t>(i, j, k, state, SYMM_INDEX(c_h11, idx1, idx2));
    }

    rhs_equation(rhs.cellData(i, j, k), state.cellData(i, j, k), d1, advec, diff2_lapse, diff2_chi, diff2_shift, diff2_h);
 
    m_deriv.add_dissipation(i, j, k, rhs.cellData(i, j, k), state, m_sigma, NUM_CCZ4_VARS);
}

// NOLINTBEGIN(readability-function-cognitive-complexity)
template <class gauge_t, class deriv_t>
template <class data_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
CCZ4RHS<gauge_t, deriv_t>::rhs_equation(
    const amrex::CellData<data_t> &rhs_cell_data, 
    const amrex::CellData<data_t const> &cell_data,
    const amrex::GpuArray<Tensor<1, data_t>, NUM_CCZ4_VARS> &d1,
    const amrex::GpuArray<data_t, NUM_CCZ4_VARS> &advec,
    const Tensor<2, data_t> &d2_lapse,
    const Tensor<2, data_t> &d2_chi,
    const Tensor<1, Tensor<2, data_t>> &d2_shift,
    const Tensor<2, Tensor<2, data_t>> &d2_h) const
{
    using namespace TensorAlgebra;

    const VarsWrapper<data_t const> vars(cell_data);
    
    auto h_UU  = compute_inverse_metric(vars);

    auto chris = compute_christoffel(d1, h_UU);

    Tensor<1, data_t> Z_over_chi;
    Tensor<1, data_t> Z; // NOLINT(readability-identifier-length)

    if (m_formulation == USE_BSSN)
    {
        FOR (i)
            Z_over_chi[i] = 0.0;
    }
    else
    {

        FOR (i)
            Z_over_chi[i] = 0.5 * (vars.Gamma(i) - chris.contracted[i]);
    }
    FOR (i)
        Z[i] = vars.chi() * Z_over_chi[i];

    auto ricci =
      CCZ4Geometry::compute_ricci_Z(vars, d1, d2_chi, d2_h, h_UU, chris, Z_over_chi);

    data_t divshift        = compute_divshift(d1);
    data_t Z_dot_d1lapse   = compute_dot_product(Z, d1[c_lapse]);
    data_t dlapse_dot_dchi = compute_dot_product(d1[c_lapse], d1[c_chi], h_UU);

    Tensor<2, data_t> covdtilde2lapse;
    Tensor<2, data_t> covd2lapse;
    FOR (k, l)
    {
        covdtilde2lapse[k][l] = d2_lapse[k][l];
        FOR (m)
        {
            covdtilde2lapse[k][l] -= chris.ULL[m][k][l] * d1[c_lapse][m];
        }
        covd2lapse[k][l] =
            vars.chi() * covdtilde2lapse[k][l] +
            0.5 * (d1[c_lapse][k] * d1[c_chi][l] + d1[c_chi][k] * d1[c_lapse][l] -
                   vars.h(k, l) * dlapse_dot_dchi);
    }

    data_t tr_covd2lapse = -(GR_SPACEDIM / 2.0) * dlapse_dot_dchi;
    FOR (i)
    {
        tr_covd2lapse -= vars.chi() * chris.contracted[i] * d1[c_lapse][i];
        FOR (j)
        {
            tr_covd2lapse += h_UU[i][j] * (vars.chi() * d2_lapse[i][j] +
                                           d1[c_lapse][i] * d1[c_chi][j]);
        }
    }

    Tensor<2, data_t> A_UU = raise_all(cell_data, c_A11, h_UU);

    // A^{ij} A_{ij}. - Note the abuse of the compute trace function.
    data_t tr_A2 = compute_trace(cell_data, c_A11, A_UU);
    rhs_cell_data[c_chi]      = advec[c_chi] +
      (2.0 / GR_SPACEDIM) * vars.chi() * (vars.lapse() * vars.K() - divshift);

    FOR (i, j)
    {
        rhs_cell_data[SYMM_INDEX(c_h11, i, j)]  = advec[SYMM_INDEX(c_h11, i, j)] - 2.0 * vars.lapse() * vars.A(i, j) -
                      (2.0 / GR_SPACEDIM) * vars.h(i, j) * divshift;
        FOR (k)
        {
            rhs_cell_data[SYMM_INDEX(c_h11, i, j)] +=
                vars.h(k, i) * d1[c_shift1 + k][j] + vars.h(k, j) * d1[c_shift1+ k][i];
        }
    }

    Tensor<2, data_t> Adot_TF;
    FOR (i, j)
    {
        Adot_TF[i][j] =
            -covd2lapse[i][j] + vars.chi() * vars.lapse() * ricci.LL[i][j];
    }
    make_trace_free(Adot_TF, vars, h_UU);

    FOR (i, j)
    {
        rhs_cell_data[SYMM_INDEX(c_A11, i, j)] = advec[SYMM_INDEX(c_A11, i, j)] + Adot_TF[i][j] +
                      vars.A(i, j) * (vars.lapse() * (vars.K() - 2 * vars.Theta()) -
                                      (2.0 / GR_SPACEDIM) * divshift);
        FOR (k)
        {
            rhs_cell_data[SYMM_INDEX(c_A11, i, j)] +=
                vars.A(k, i) * d1[c_shift1 + k][j] + vars.A(k, j) * d1[c_shift1 + k][i];
            FOR (l)
            {
                rhs_cell_data[SYMM_INDEX(c_A11, i, j)] -=
                    2 * vars.lapse() * h_UU[k][l] * vars.A(i, k) * vars.A(l, j);
            }
        }
    }

    data_t kappa1_times_lapse;
    if (m_params.covariantZ4)
    {
        kappa1_times_lapse = m_params.kappa1;
    }
    else
    {
        kappa1_times_lapse = m_params.kappa1 * vars.lapse();
    }

    if (m_formulation == USE_BSSN)
    {
        rhs_cell_data[c_Theta] = 0; // ensure the Theta of CCZ4 remains at zero
        // Use hamiltonian constraint to remove ricci.scalar for BSSN update
        rhs_cell_data[c_K] = advec[c_K] + vars.lapse() * (tr_A2 + vars.K() * vars.K() / GR_SPACEDIM) -
                tr_covd2lapse;
        rhs_cell_data[c_K] += -2 * vars.lapse() * m_cosmological_constant / (GR_SPACEDIM - 1.);
    }
    else
    {
        rhs_cell_data[c_Theta] =
            advec[c_Theta] +
            0.5 * vars.lapse() *
                (ricci.scalar - tr_A2 +
                 ((GR_SPACEDIM - 1.0) / (double)GR_SPACEDIM) * vars.K() * vars.K() -
                 2 * vars.Theta() * vars.K()) -
            0.5 * vars.Theta() * kappa1_times_lapse *
                ((GR_SPACEDIM + 1) + m_params.kappa2 * (GR_SPACEDIM - 1)) -
            Z_dot_d1lapse;

        rhs_cell_data[c_Theta] += -vars.lapse() * m_cosmological_constant;
        rhs_cell_data[c_K] =
            advec[c_K] +
            vars.lapse() * (ricci.scalar + vars.K() * (vars.K() - 2 * vars.Theta())) -
            kappa1_times_lapse * GR_SPACEDIM * (1 + m_params.kappa2) *
                vars.Theta() -
            tr_covd2lapse;
        rhs_cell_data[c_K] += -2 * vars.lapse() * GR_SPACEDIM / (GR_SPACEDIM - 1.) *
                 m_cosmological_constant;
    }

    Tensor<1, data_t> Gammadot;
    FOR (i)
    {
        Gammadot[i] = (2.0 / GR_SPACEDIM) *
                          (divshift * (chris.contracted[i] +
                                       2 * m_params.kappa3 * Z_over_chi[i]) -
                           2 * vars.lapse() * vars.K() * Z_over_chi[i]) -
                      2 * kappa1_times_lapse * Z_over_chi[i];
        FOR (j)
        {
            Gammadot[i] +=
                2 * h_UU[i][j] *
                    (vars.lapse() * d1[c_Theta][j] - vars.Theta() * d1[c_lapse][j]) -
                2 * A_UU[i][j] * d1[c_lapse][j] -
                vars.lapse() * ((2 * (GR_SPACEDIM - 1.0) / (double)GR_SPACEDIM) *
                                  h_UU[i][j] * d1[c_K][j] +
                              GR_SPACEDIM * A_UU[i][j] * d1[c_chi][j] / vars.chi()) -
                (chris.contracted[j] + 2 * m_params.kappa3 * Z_over_chi[j]) *
                    d1[c_shift1+ i][j];

            FOR (k)
            {
                Gammadot[i] +=
                    2 * vars.lapse() * chris.ULL[i][j][k] * A_UU[j][k] +
                    h_UU[j][k] * d2_shift[i][j][k] +
                    ((GR_SPACEDIM - 2.0) / (double)GR_SPACEDIM) * h_UU[i][j] *
                        d2_shift[k][j][k];
            }
        }
    }

    Tensor<1, data_t> rhs_Gamma;
    FOR (i)
    {
        rhs_cell_data[c_Gamma1 + i] = advec[c_Gamma1 + i] + Gammadot[i];
    }

    m_gauge.rhs_gauge(rhs_cell_data, cell_data, d1, advec);
}
// NOLINTEND(readability-function-cognitive-complexity)

#endif /* CCZ4RHS_IMPL_HPP_ */
