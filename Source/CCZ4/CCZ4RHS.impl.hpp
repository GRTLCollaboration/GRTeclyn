/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(CCZ4RHS_HPP_)
#error "This file should only be included through CCZ4RHS.hpp"
#endif

#ifndef CCZ4RHS_IMPL_HPP_
#define CCZ4RHS_IMPL_HPP_

#include "DimensionDefinitions.hpp"
#include "GRInterval.hpp"
#include "VarsTools.hpp"

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
    const auto vars = load_vars<Vars>(state.cellData(i, j, k));
    const auto d1   = m_deriv.template diff1<Vars>(i, j, k, state);
    const auto d2   = m_deriv.template diff2<Diff2Vars>(i, j, k, state);
    const auto advec =
        m_deriv.template advection<Vars>(i, j, k, state, vars.shift);

    Vars<data_t> rhs_vars;
    rhs_equation(rhs_vars, state.cellData(i, j, k), d1, d2, advec);

    m_deriv.add_dissipation(i, j, k, rhs_vars, state, m_sigma);

    store_vars(rhs.cellData(i, j, k), rhs_vars);
}

// NOLINTBEGIN(readability-function-cognitive-complexity)
template <class gauge_t, class deriv_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
CCZ4RHS<gauge_t, deriv_t>::rhs_equation(
    vars_t<data_t> &rhs, const amrex::CellData<data_t const> &cell_data,
    const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2,
    const vars_t<data_t> &advec) const
{
    using namespace TensorAlgebra;

    Tensor<2,data_t> h; //{cell_data[c_h11], cell_data[c_h12], cell_data[c_h13], cell_data[c_h22], cell_data[c_h23], cell_data[c_h33]};
    h[0][0] = cell_data[c_h11]; // change to 1D
    h[0][1] = cell_data[c_h12];
    h[0][2] = cell_data[c_h13];
    h[1][0] = cell_data[c_h12];
    h[1][1] = cell_data[c_h22];
    h[1][2] = cell_data[c_h23];
    h[2][0] = cell_data[c_h13];
    h[2][1] = cell_data[c_h23];
    h[2][2] = cell_data[c_h33];


    Tensor<1, data_t> Gamma;
    Gamma[0] = cell_data[c_Gamma1];
    Gamma[1] = cell_data[c_Gamma2];
    Gamma[2] = cell_data[c_Gamma3];

    Tensor<1, data_t> shift;
    shift[0] = cell_data[c_shift1];
    shift[1] = cell_data[c_shift2];
    shift[2] = cell_data[c_shift3];

    Tensor<2,data_t> A;
    A[0][0] = cell_data[c_A11]; // change to 1D
    A[0][1] = cell_data[c_A12];
    A[0][2] = cell_data[c_A13];
    A[1][0] = cell_data[c_A12];
    A[1][1] = cell_data[c_A22];
    A[1][2] = cell_data[c_A23];
    A[2][0] = cell_data[c_A13];
    A[2][1] = cell_data[c_A23];
    A[2][2] = cell_data[c_A33];
    
    auto h_UU  = compute_inverse_sym(h);
    auto chris = compute_christoffel(d1.h, h_UU);

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
            Z_over_chi[i] = 0.5 * (Gamma[i] - chris.contracted[i]);
    }
    FOR (i)
        Z[i] = cell_data[c_chi] * Z_over_chi[i];

    auto ricci =
      CCZ4Geometry::compute_ricci_Z(cell_data[c_chi], h, Gamma, d1, d2, h_UU, chris, Z_over_chi);

    data_t divshift        = compute_trace(d1.shift);
    data_t Z_dot_d1lapse   = compute_dot_product(Z, d1.lapse);
    data_t dlapse_dot_dchi = compute_dot_product(d1.lapse, d1.chi, h_UU);

    Tensor<2, data_t> covdtilde2lapse;
    Tensor<2, data_t> covd2lapse;
    FOR (k, l)
    {
        covdtilde2lapse[k][l] = d2.lapse[k][l];
        FOR (m)
        {
            covdtilde2lapse[k][l] -= chris.ULL[m][k][l] * d1.lapse[m];
        }
        covd2lapse[k][l] =
            cell_data[c_chi] * covdtilde2lapse[k][l] +
            0.5 * (d1.lapse[k] * d1.chi[l] + d1.chi[k] * d1.lapse[l] -
                   h[k][l] * dlapse_dot_dchi);
    }

    data_t tr_covd2lapse = -(GR_SPACEDIM / 2.0) * dlapse_dot_dchi;
    FOR (i)
    {
        tr_covd2lapse -= cell_data[c_chi] * chris.contracted[i] * d1.lapse[i];
        FOR (j)
        {
            tr_covd2lapse += h_UU[i][j] * (cell_data[c_chi] * d2.lapse[i][j] +
                                           d1.lapse[i] * d1.chi[j]);
        }
    }

    Tensor<2, data_t> A_UU = raise_all(A, h_UU);

    // A^{ij} A_{ij}. - Note the abuse of the compute trace function.
    data_t tr_A2 = compute_trace(A, A_UU);
    rhs.chi      = advec.chi +
      (2.0 / GR_SPACEDIM) * cell_data[c_chi] * (cell_data[c_lapse] * cell_data[c_K] - divshift);
    FOR (i, j)
    {
        rhs.h[i][j] = advec.h[i][j] - 2.0 * cell_data[c_lapse] * A[i][j] -
                      (2.0 / GR_SPACEDIM) * h[i][j] * divshift;
        FOR (k)
        {
            rhs.h[i][j] +=
                h[k][i] * d1.shift[k][j] + h[k][j] * d1.shift[k][i];
        }
    }

    Tensor<2, data_t> Adot_TF;
    FOR (i, j)
    {
        Adot_TF[i][j] =
            -covd2lapse[i][j] + cell_data[c_chi] * cell_data[c_lapse] * ricci.LL[i][j];
    }
    make_trace_free(Adot_TF, h, h_UU);

    FOR (i, j)
    {
        rhs.A[i][j] = advec.A[i][j] + Adot_TF[i][j] +
                      A[i][j] * (cell_data[c_lapse] * (cell_data[c_K] - 2 * cell_data[c_Theta]) -
                                      (2.0 / GR_SPACEDIM) * divshift);
        FOR (k)
        {
            rhs.A[i][j] +=
                A[k][i] * d1.shift[k][j] + A[k][j] * d1.shift[k][i];
            FOR (l)
            {
                rhs.A[i][j] -=
                    2 * cell_data[c_lapse] * h_UU[k][l] * A[i][k] * A[l][j];
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
        kappa1_times_lapse = m_params.kappa1 * cell_data[c_lapse];
    }

    if (m_formulation == USE_BSSN)
    {
        rhs.Theta = 0; // ensure the Theta of CCZ4 remains at zero
        // Use hamiltonian constraint to remove ricci.scalar for BSSN update
        rhs.K = advec.K + cell_data[c_lapse] * (tr_A2 + cell_data[c_K] * cell_data[c_K] / GR_SPACEDIM) -
                tr_covd2lapse;
        rhs.K += -2 * cell_data[c_lapse] * m_cosmological_constant / (GR_SPACEDIM - 1.);
    }
    else
    {
        rhs.Theta =
            advec.Theta +
            0.5 * cell_data[c_lapse] *
                (ricci.scalar - tr_A2 +
                 ((GR_SPACEDIM - 1.0) / (double)GR_SPACEDIM) * cell_data[c_K] * cell_data[c_K] -
                 2 * cell_data[c_Theta] * cell_data[c_K]) -
            0.5 * cell_data[c_Theta] * kappa1_times_lapse *
                ((GR_SPACEDIM + 1) + m_params.kappa2 * (GR_SPACEDIM - 1)) -
            Z_dot_d1lapse;

        rhs.Theta += -cell_data[c_lapse] * m_cosmological_constant;
        rhs.K =
            advec.K +
            cell_data[c_lapse] * (ricci.scalar + cell_data[c_K] * (cell_data[c_K] - 2 * cell_data[c_Theta])) -
            kappa1_times_lapse * GR_SPACEDIM * (1 + m_params.kappa2) *
                cell_data[c_Theta] -
            tr_covd2lapse;
        rhs.K += -2 * cell_data[c_lapse] * GR_SPACEDIM / (GR_SPACEDIM - 1.) *
                 m_cosmological_constant;
    }

    Tensor<1, data_t> Gammadot;
    FOR (i)
    {
        Gammadot[i] = (2.0 / GR_SPACEDIM) *
                          (divshift * (chris.contracted[i] +
                                       2 * m_params.kappa3 * Z_over_chi[i]) -
                           2 * cell_data[c_lapse] * cell_data[c_K] * Z_over_chi[i]) -
                      2 * kappa1_times_lapse * Z_over_chi[i];
        FOR (j)
        {
            Gammadot[i] +=
                2 * h_UU[i][j] *
                    (cell_data[c_lapse] * d1.Theta[j] - cell_data[c_Theta] * d1.lapse[j]) -
                2 * A_UU[i][j] * d1.lapse[j] -
                cell_data[c_lapse] * ((2 * (GR_SPACEDIM - 1.0) / (double)GR_SPACEDIM) *
                                  h_UU[i][j] * d1.K[j] +
                              GR_SPACEDIM * A_UU[i][j] * d1.chi[j] / cell_data[c_chi]) -
                (chris.contracted[j] + 2 * m_params.kappa3 * Z_over_chi[j]) *
                    d1.shift[i][j];

            FOR (k)
            {
                Gammadot[i] +=
                    2 * cell_data[c_lapse] * chris.ULL[i][j][k] * A_UU[j][k] +
                    h_UU[j][k] * d2.shift[i][j][k] +
                    ((GR_SPACEDIM - 2.0) / (double)GR_SPACEDIM) * h_UU[i][j] *
                        d2.shift[k][j][k];
            }
        }
    }

    FOR (i)
    {
        rhs.Gamma[i] = advec.Gamma[i] + Gammadot[i];
    }


    Tensor<1,data_t> B;
    B[0] = cell_data[c_B1];
    B[1] = cell_data[c_B2];
    B[2] = cell_data[c_B3]; 
    //need K, lapse, B, Theta
    m_gauge.rhs_gauge(rhs, cell_data[c_lapse], cell_data[c_K], cell_data[c_Theta], B, d1, d2, advec);
}
// NOLINTEND(readability-function-cognitive-complexity)

#endif /* CCZ4RHS_IMPL_HPP_ */
