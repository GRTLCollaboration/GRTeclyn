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
    Tensor<1, data_t> shift;
    shift[0] = state.cellData(i, j, k)[c_shift1];
    shift[1] = state.cellData(i, j, k)[c_shift2];
    shift[2] = state.cellData(i, j, k)[c_shift3];

    const auto advec =
        m_deriv.template advection<data_t, NUM_CCZ4_VARS>(i, j, k, state, shift);

    const auto d1   = m_deriv.template diff1<data_t, NUM_CCZ4_VARS>(i, j, k, state);

    Tensor<2, data_t> diff2_lapse = m_deriv.template diff2<data_t>(i, j, k, state, c_lapse);
    Tensor<2, data_t> diff2_chi = m_deriv.template diff2<data_t>(i, j, k, state, c_chi);

    Tensor<1, Tensor<2, data_t>> diff2_shift;
    diff2_shift[0] = m_deriv.template diff2<data_t>(i, j, k, state, c_shift1);
    diff2_shift[1] = m_deriv.template diff2<data_t>(i, j, k, state, c_shift2);
    diff2_shift[2] = m_deriv.template diff2<data_t>(i, j, k, state, c_shift3);

    Tensor<2, Tensor<2, data_t>> diff2_h;
    diff2_h[0][0] = m_deriv.template diff2<data_t>(i, j, k, state, c_h11);
    diff2_h[0][1] = m_deriv.template diff2<data_t>(i, j, k, state, c_h12);
    diff2_h[0][2] = m_deriv.template diff2<data_t>(i, j, k, state, c_h13);
    diff2_h[1][0] = m_deriv.template diff2<data_t>(i, j, k, state, c_h12);
    diff2_h[1][1] = m_deriv.template diff2<data_t>(i, j, k, state, c_h22);
    diff2_h[1][2] = m_deriv.template diff2<data_t>(i, j, k, state, c_h23);
    diff2_h[2][0] = m_deriv.template diff2<data_t>(i, j, k, state, c_h13);
    diff2_h[2][1] = m_deriv.template diff2<data_t>(i, j, k, state, c_h23);
    diff2_h[2][2] = m_deriv.template diff2<data_t>(i, j, k, state, c_h33);

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

    Tensor<2,data_t> h;
    h[0][0] = cell_data[c_h11]; // change to 1D
    h[0][1] = cell_data[c_h12];
    h[0][2] = cell_data[c_h13];
    h[1][0] = cell_data[c_h12];
    h[1][1] = cell_data[c_h22];
    h[1][2] = cell_data[c_h23];
    h[2][0] = cell_data[c_h13];
    h[2][1] = cell_data[c_h23];
    h[2][2] = cell_data[c_h33];
    
    auto h_UU  = compute_inverse_sym(h);
    // FIXME: d1.h is no longer immediately a Tensor object

    Tensor<2, Tensor<1,data_t>> d1_h; //{cell_data[c_h11], cell_data[c_h12], cell_data[c_h13], cell_data[c_h22], cell_data[c_h23], cell_data[c_h33]};
    Tensor<2,data_t> d1_shift;
    FOR (i)
    {
        d1_h[0][0][i] = d1[c_h11][i]; // change to 1D
        d1_h[0][1][i] = d1[c_h12][i];
        d1_h[0][2][i] = d1[c_h13][i];
        d1_h[1][0][i] = d1[c_h12][i];
        d1_h[1][1][i] = d1[c_h22][i];
        d1_h[1][2][i] = d1[c_h23][i];
        d1_h[2][0][i] = d1[c_h13][i];
        d1_h[2][1][i] = d1[c_h23][i];
        d1_h[2][2][i] = d1[c_h33][i];

        d1_shift[0][i] = d1[c_shift1][i];
        d1_shift[1][i] = d1[c_shift2][i];
        d1_shift[2][i] = d1[c_shift3][i];
    }
    

    auto chris = compute_christoffel(d1_h, h_UU);

    Tensor<1, data_t> Z_over_chi;
    Tensor<1, data_t> Z; // NOLINT(readability-identifier-length)

    Tensor<1, data_t> Gamma;
    Gamma[0] = cell_data[c_Gamma1];
    Gamma[1] = cell_data[c_Gamma2];
    Gamma[2] = cell_data[c_Gamma3];

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
      CCZ4Geometry::compute_ricci_Z(cell_data[c_chi], h, Gamma, d1, d2_chi, d2_h, h_UU, chris, Z_over_chi, d1_h, true);

    data_t divshift        = compute_trace(d1_shift);
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
            cell_data[c_chi] * covdtilde2lapse[k][l] +
            0.5 * (d1[c_lapse][k] * d1[c_chi][l] + d1[c_chi][k] * d1[c_lapse][l] -
                   h[k][l] * dlapse_dot_dchi);
    }

    data_t tr_covd2lapse = -(GR_SPACEDIM / 2.0) * dlapse_dot_dchi;
    FOR (i)
    {
        tr_covd2lapse -= cell_data[c_chi] * chris.contracted[i] * d1[c_lapse][i];
        FOR (j)
        {
            tr_covd2lapse += h_UU[i][j] * (cell_data[c_chi] * d2_lapse[i][j] +
                                           d1[c_lapse][i] * d1[c_chi][j]);
        }
    }

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

    Tensor<2, data_t> A_UU = raise_all(A, h_UU);

    // A^{ij} A_{ij}. - Note the abuse of the compute trace function.
    data_t tr_A2 = compute_trace(A, A_UU);
    rhs_cell_data[c_chi]      = advec[c_chi] +
      (2.0 / GR_SPACEDIM) * cell_data[c_chi] * (cell_data[c_lapse] * cell_data[c_K] - divshift);

    Tensor<2,data_t> rhs_h;
    /*Tensor<2,data_t> advec_h;

    advec_h[0][0] = advec[c_h11]; // change to 1D
    advec_h[0][1] = advec[c_h12];
    advec_h[0][2] = advec[c_h13];
    advec_h[1][0] = advec[c_h12];
    advec_h[1][1] = advec[c_h22];
    advec_h[1][2] = advec[c_h23];
    advec_h[2][0] = advec[c_h13];
    advec_h[2][1] = advec[c_h23];
    advec_h[2][2] = advec[c_h33];*/

    FOR (i, j)
    {
        rhs_h[i][j] = advec[SYMM_INDEX(c_h11, i, j)] - 2.0 * cell_data[c_lapse] * A[i][j] -
                      (2.0 / GR_SPACEDIM) * h[i][j] * divshift;
        FOR (k)
        {
            rhs_h[i][j] +=
                h[k][i] * d1_shift[k][j] + h[k][j] * d1_shift[k][i];
        }
    }

    Tensor<2, data_t> Adot_TF;
    FOR (i, j)
    {
        Adot_TF[i][j] =
            -covd2lapse[i][j] + cell_data[c_chi] * cell_data[c_lapse] * ricci.LL[i][j];
    }
    make_trace_free(Adot_TF, h, h_UU);

    Tensor<2,data_t> rhs_A;
    /*Tensor<2,data_t> advec_A;

    advec_A[0][0] = advec[c_A11]; // change to 1D
    advec_A[0][1] = advec[c_A12];
    advec_A[0][2] = advec[c_A13];
    advec_A[1][0] = advec[c_A12];
    advec_A[1][1] = advec[c_A22];
    advec_A[1][2] = advec[c_A23];
    advec_A[2][0] = advec[c_A13];
    advec_A[2][1] = advec[c_A23];
    advec_A[2][2] = advec[c_A33];*/
    FOR (i, j)
    {
        rhs_A[i][j] = advec[SYMM_INDEX(c_A11, i, j)] + Adot_TF[i][j] +
                      A[i][j] * (cell_data[c_lapse] * (cell_data[c_K] - 2 * cell_data[c_Theta]) -
                                      (2.0 / GR_SPACEDIM) * divshift);
        FOR (k)
        {
            rhs_A[i][j] +=
                A[k][i] * d1_shift[k][j] + A[k][j] * d1_shift[k][i];
            FOR (l)
            {
                rhs_A[i][j] -=
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
        rhs_cell_data[c_Theta] = 0; // ensure the Theta of CCZ4 remains at zero
        // Use hamiltonian constraint to remove ricci.scalar for BSSN update
        rhs_cell_data[c_K] = advec[c_K] + cell_data[c_lapse] * (tr_A2 + cell_data[c_K] * cell_data[c_K] / GR_SPACEDIM) -
                tr_covd2lapse;
        rhs_cell_data[c_K] += -2 * cell_data[c_lapse] * m_cosmological_constant / (GR_SPACEDIM - 1.);
    }
    else
    {
        rhs_cell_data[c_Theta] =
            advec[c_Theta] +
            0.5 * cell_data[c_lapse] *
                (ricci.scalar - tr_A2 +
                 ((GR_SPACEDIM - 1.0) / (double)GR_SPACEDIM) * cell_data[c_K] * cell_data[c_K] -
                 2 * cell_data[c_Theta] * cell_data[c_K]) -
            0.5 * cell_data[c_Theta] * kappa1_times_lapse *
                ((GR_SPACEDIM + 1) + m_params.kappa2 * (GR_SPACEDIM - 1)) -
            Z_dot_d1lapse;

        rhs_cell_data[c_Theta] += -cell_data[c_lapse] * m_cosmological_constant;
        rhs_cell_data[c_K] =
            advec[c_K] +
            cell_data[c_lapse] * (ricci.scalar + cell_data[c_K] * (cell_data[c_K] - 2 * cell_data[c_Theta])) -
            kappa1_times_lapse * GR_SPACEDIM * (1 + m_params.kappa2) *
                cell_data[c_Theta] -
            tr_covd2lapse;
        rhs_cell_data[c_K] += -2 * cell_data[c_lapse] * GR_SPACEDIM / (GR_SPACEDIM - 1.) *
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
                    (cell_data[c_lapse] * d1[c_Theta][j] - cell_data[c_Theta] * d1[c_lapse][j]) -
                2 * A_UU[i][j] * d1[c_lapse][j] -
                cell_data[c_lapse] * ((2 * (GR_SPACEDIM - 1.0) / (double)GR_SPACEDIM) *
                                  h_UU[i][j] * d1[c_K][j] +
                              GR_SPACEDIM * A_UU[i][j] * d1[c_chi][j] / cell_data[c_chi]) -
                (chris.contracted[j] + 2 * m_params.kappa3 * Z_over_chi[j]) *
                    d1_shift[i][j];

            FOR (k)
            {
                Gammadot[i] +=
                    2 * cell_data[c_lapse] * chris.ULL[i][j][k] * A_UU[j][k] +
                    h_UU[j][k] * d2_shift[i][j][k] +
                    ((GR_SPACEDIM - 2.0) / (double)GR_SPACEDIM) * h_UU[i][j] *
                        d2_shift[k][j][k];
            }
        }
    }

    Tensor<1, data_t> rhs_Gamma;
    FOR (i)
    {
        rhs_Gamma[i] = advec[c_Gamma1 + i] + Gammadot[i];
    }

    rhs_cell_data[c_Gamma1]= rhs_Gamma[0];
    rhs_cell_data[c_Gamma2]= rhs_Gamma[1];
    rhs_cell_data[c_Gamma3]= rhs_Gamma[2];

    rhs_cell_data[c_h11]=rhs_h[0][0]; // change to 1D
    rhs_cell_data[c_h12]=rhs_h[0][1];
    rhs_cell_data[c_h13]=rhs_h[0][2];
    rhs_cell_data[c_h22]=rhs_h[1][1];
    rhs_cell_data[c_h23]=rhs_h[1][2];
    rhs_cell_data[c_h33]=rhs_h[2][2];

    rhs_cell_data[c_A11]=rhs_A[0][0]; // change to 1D
    rhs_cell_data[c_A12]=rhs_A[0][1];
    rhs_cell_data[c_A13]=rhs_A[0][2];
    rhs_cell_data[c_A22]=rhs_A[1][1];
    rhs_cell_data[c_A23]=rhs_A[2][1];
    rhs_cell_data[c_A33]=rhs_A[2][2];


    Tensor<1,data_t> B;
    B[0] = cell_data[c_B1];
    B[1] = cell_data[c_B2];
    B[2] = cell_data[c_B3]; 

    //need K, lapse, B, Theta
    m_gauge.rhs_gauge(rhs_cell_data, cell_data[c_lapse], cell_data[c_K], cell_data[c_Theta], B, d1, advec);
}
// NOLINTEND(readability-function-cognitive-complexity)

#endif /* CCZ4RHS_IMPL_HPP_ */
