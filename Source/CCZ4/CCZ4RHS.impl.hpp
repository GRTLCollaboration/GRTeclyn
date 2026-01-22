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
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CCZ4RHS<gauge_t, deriv_t>::operator()(
    int ix, int iy, int iz, const amrex::Array4<amrex::Real> &rhs_state,
    const amrex::Array4<amrex::Real const> &state) const
{
    const amrex::CellData<const amrex::Real> &state_cell_data =
        state.cellData(ix, iy, iz);
    CCZ4Vars vars(state_cell_data);

    //    Get the derivatives
    const CCZ4D1Vars d1(ix, iy, iz, state, m_deriv);
    const CCZ4D2Vars d2(ix, iy, iz, state, m_deriv);
    const CCZ4AdvecVars advec(ix, iy, iz, state, m_deriv);

    const amrex::CellData<amrex::Real> &rhs_cell_data =
        rhs_state.cellData(ix, iy, iz);
    // CCZ4Vars rhs(rhs_cell_data);

    rhs_equation(rhs_cell_data, vars, d1, d2, advec);

    // calculate_chi_rhs(ix, iy, iz, rhs_state, state);

    // calculate_h_ij_rhs(ix, iy, iz, rhs_state, state);

    // calculate_A_ij_rhs(ix, iy, iz, rhs_state, state);

    // calculate_Theta_rhs(ix, iy, iz, rhs_state, state);

    // calculate_Gamma_rhs(ix, iy, iz, rhs_state, state);

    // apply_gauge(ix, iy, iz, rhs_state, state);

    m_deriv.add_dissipation(ix, iy, iz, rhs_cell_data, state, m_sigma,
                            NUM_CCZ4_VARS);
}

template <class gauge_t, class deriv_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
CCZ4RHS<gauge_t, deriv_t>::calculate_chi_rhs(
    int ix, int iy, int iz, const amrex::Array4<amrex::Real> &rhs_state,
    const amrex::Array4<amrex::Real const> &state) const
{
    const amrex::CellData<amrex::Real> &rhs_cell_data =
        rhs_state.cellData(ix, iy, iz);
    const amrex::CellData<const amrex::Real> &state_cell_data =
        state.cellData(ix, iy, iz);

    Tensor<1, amrex::Real> shift_vector;
    FOR (idir)
    {
        shift_vector[idir] = state(ix, iy, iz, c_shift1 + idir);
    }

    Tensor<2, amrex::Real> d1_shift =
        m_deriv.diff1_vector(ix, iy, iz, state, c_shift1);
    amrex::Real divshift = CCZ4Geometry::compute_divshift(d1_shift);
    amrex::Real advec_chi =
        m_deriv.advection(ix, iy, iz, state, shift_vector, c_chi);
    rhs_cell_data[c_chi] =
        advec_chi +
        (2.0 / (double)GR_SPACEDIM) * state_cell_data[c_chi] *
            (state_cell_data[c_lapse] * state_cell_data[c_K] - divshift);

    // Calculation of h_ij RHS

    Tensor<2, amrex::Real> advec_h =
        m_deriv.advec_tensor(ix, iy, iz, state, shift_vector, c_h11);

    FOR (i, j)
    {
        rhs_cell_data[var_idx(c_h11, i, j)] =
            advec_h[i][j] -
            2.0 * state_cell_data[c_lapse] *
                state_cell_data[var_idx(c_A11, i, j)] -
            (2.0 / (double)GR_SPACEDIM) *
                state_cell_data[var_idx(c_h11, i, j)] * divshift;

        FOR (k)
        {
            rhs_cell_data[var_idx(c_h11, i, j)] +=
                state_cell_data[var_idx(c_h11, k, i)] * d1_shift[k][j] +
                state_cell_data[var_idx(c_h11, k, j)] * d1_shift[k][i];
        }
    }
}

template <class gauge_t, class deriv_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
CCZ4RHS<gauge_t, deriv_t>::calculate_h_ij_rhs(
    int ix, int iy, int iz, const amrex::Array4<amrex::Real> &rhs_state,
    const amrex::Array4<amrex::Real const> &state) const
{

    const amrex::CellData<amrex::Real> &rhs_cell_data =
        rhs_state.cellData(ix, iy, iz);
    const amrex::CellData<const amrex::Real> &state_cell_data =
        state.cellData(ix, iy, iz);

    Tensor<1, amrex::Real> shift_vector;
    FOR (idir)
    {
        shift_vector[idir] = state(ix, iy, iz, c_shift1 + idir);
    }

    Tensor<2, amrex::Real> d1_shift =
        m_deriv.diff1_vector(ix, iy, iz, state, c_shift1);
    amrex::Real divshift = CCZ4Geometry::compute_divshift(d1_shift);
    Tensor<2, amrex::Real> advec_h =
        m_deriv.advec_tensor(ix, iy, iz, state, shift_vector, c_h11);

    FOR (i, j)
    {
        rhs_cell_data[var_idx(c_h11, i, j)] =
            advec_h[i][j] -
            2.0 * state_cell_data[c_lapse] *
                state_cell_data[var_idx(c_A11, i, j)] -
            (2.0 / (double)GR_SPACEDIM) *
                state_cell_data[var_idx(c_h11, i, j)] * divshift;

        FOR (k)
        {
            rhs_cell_data[var_idx(c_h11, i, j)] +=
                state_cell_data[var_idx(c_h11, k, i)] * d1_shift[k][j] +
                state_cell_data[var_idx(c_h11, k, j)] * d1_shift[k][i];
        }
    }
}

template <class gauge_t, class deriv_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
CCZ4RHS<gauge_t, deriv_t>::calculate_A_ij_rhs_use_amrex_array(
    int ix, int iy, int iz, const amrex::Array4<amrex::Real> &rhs_state,
    const amrex::Array4<amrex::Real const> &state) const
{
    const amrex::CellData<amrex::Real> &rhs_cell_data =
        rhs_state.cellData(ix, iy, iz);
    const amrex::CellData<const amrex::Real> &state_cell_data =
        state.cellData(ix, iy, iz);

    const auto h_UU =
        CCZ4Geometry::compute_inverse_metric_array(state_cell_data);

    // hij derivatives
    auto d1_h = m_deriv.diff1_array_tensor(ix, iy, iz, state, c_h11);

    const auto chris = CCZ4Geometry::compute_christoffel(d1_h, h_UU);

    Tensor<1, amrex::Real> Z_over_chi;
    Tensor<1, amrex::Real> Z; // NOLINT(readability-identifier-length)

    if (m_formulation == USE_BSSN)
    {
        FOR (i)
            Z_over_chi[i] = 0.0;
    }
    else
    {

        FOR (i)
            Z_over_chi[i] =
                0.5 * (state_cell_data[c_Gamma1 + i] - chris.contracted(i));
    }
    FOR (i)
        Z[i] = state_cell_data[c_chi] * Z_over_chi[i];

    auto ricci = CCZ4Geometry::compute_ricci_Z(
        ix, iy, iz, rhs_state, state, h_UU, chris, Z_over_chi, m_deriv);

    Tensor<2, amrex::Real> d1_shift =
        m_deriv.diff1_vector(ix, iy, iz, state, c_shift1);
    amrex::Real divshift = CCZ4Geometry::compute_divshift(d1_shift);

    auto d1_lapse = m_deriv.diff1_scalar(ix, iy, iz, state, c_lapse);
    auto d1_chi   = m_deriv.diff1_scalar(ix, iy, iz, state, c_chi);

    amrex::Real Z_dot_d1lapse = TensorAlgebra::compute_dot_product(Z, d1_lapse);
    amrex::Real dlapse_dot_dchi =
        TensorAlgebra::compute_dot_product(d1_lapse, d1_chi, h_UU);

    Tensor<2, amrex::Real> covdtilde2lapse;
    Tensor<2, amrex::Real> covd2lapse;
    auto d2_lapse = m_deriv.diff2_sym_scalar(ix, iy, iz, state, c_lapse);

    FOR (k, l)
    {
        int idx               = k + l + ((k * l != 0) ? 1 : 0);
        covdtilde2lapse[k][l] = d2_lapse(idx);
        FOR (m)
        {
            covdtilde2lapse[k][l] -= chris.ULL(m, k, l) * d1_lapse[m];
        }
        covd2lapse[k][l] =
            state_cell_data[c_chi] * covdtilde2lapse[k][l] +
            0.5 * (d1_lapse[k] * d1_chi[l] + d1_chi[k] * d1_lapse[l] -
                   state_cell_data[var_idx(c_h11, k, l)] * dlapse_dot_dchi);
    }

    auto A_UU = CCZ4Geometry::compute_A_UU(state_cell_data, h_UU);

    amrex::Array2D<amrex::Real, 0, 3, 0, 3> Adot_TF;
    FOR (i, j)
    {
        Adot_TF(i, j) = -covd2lapse[i][j] + state_cell_data[c_chi] *
                                                state_cell_data[c_lapse] *
                                                ricci.LL(i, j);
    }
    CCZ4Geometry::make_trace_free(Adot_TF, state_cell_data, h_UU);

    Tensor<1, amrex::Real> shift_vector;
    FOR (idir)
    {
        shift_vector[idir] = state(ix, iy, iz, c_shift1 + idir);
    }
    auto advec_A = m_deriv.advec_tensor(ix, iy, iz, state, shift_vector, c_A11);

    FOR (i, j)
    {
        int idx = var_idx(c_A11, i, j);
        rhs_cell_data[idx] =
            advec_A[i][j] + Adot_TF(i, j) +
            state_cell_data[idx] *
                (state_cell_data[c_lapse] *
                     (state_cell_data[c_K] - 2.0 * state_cell_data[c_Theta]) -
                 (2.0 / (double)GR_SPACEDIM) * divshift);
        FOR (k)
        {
            rhs_cell_data[idx] +=
                state_cell_data[var_idx(c_A11, k, i)] * d1_shift[k][j] +
                state_cell_data[var_idx(c_A11, k, j)] * d1_shift[k][i];
            FOR (l)
            {
                rhs_cell_data[idx] -= 2.0 * state_cell_data[c_lapse] *
                                      h_UU(k, l) *
                                      state_cell_data[var_idx(c_A11, i, k)] *
                                      state_cell_data[var_idx(c_A11, l, j)];
            }
        }
    }

    // Theta specific parts

    amrex::Real Aij_squared =
        CCZ4Geometry::compute_Aij_squared(state_cell_data, A_UU);

    amrex::Real advec_K =
        m_deriv.advection(ix, iy, iz, state, shift_vector, c_K);

    amrex::Real tr_covd2lapse = -((double)GR_SPACEDIM / 2.0) * dlapse_dot_dchi;
    FOR (i)
    {
        tr_covd2lapse -=
            state_cell_data[c_chi] * chris.contracted(i) * d1_lapse[i];
        FOR (j)
        {
            tr_covd2lapse += h_UU(i, j) * (state_cell_data[c_chi] *
                                               d2_lapse(SYMM_IDX(i, j)) +
                                           d1_lapse[i] * d1_chi[j]);
        }
    }

    amrex::Real kappa1_times_lapse;
    if (m_params.covariantZ4)
    {
        kappa1_times_lapse = m_params.kappa1;
    }
    else
    {
        kappa1_times_lapse = m_params.kappa1 * state_cell_data[c_lapse];
    }

    if (m_formulation == USE_BSSN)
    {
        // ensure the Theta of CCZ4 remains at zero
        rhs_cell_data[c_Theta] = 0.0;
        // Use hamiltonian constraint to remove ricci.scalar for BSSN update
        rhs_cell_data[c_K] =
            advec_K +
            state_cell_data[c_lapse] *
                (Aij_squared + state_cell_data[c_K] * state_cell_data[c_K] /
                                   (double)GR_SPACEDIM) -
            tr_covd2lapse -
            2.0 * state_cell_data[c_lapse] * m_cosmological_constant /
                ((double)GR_SPACEDIM - 1.0);
    }
    else
    {
        amrex::Real advec_Theta =
            m_deriv.advection(ix, iy, iz, state, shift_vector, c_Theta);

        rhs_cell_data[c_Theta] =
            advec_Theta +
            0.5 * state_cell_data[c_lapse] *
                (ricci.scalar - Aij_squared +
                 (((double)GR_SPACEDIM - 1.0) / (double)GR_SPACEDIM) *
                     state_cell_data[c_K] * state_cell_data[c_K] -
                 2.0 * state_cell_data[c_Theta] * state_cell_data[c_K]) -
            0.5 * state_cell_data[c_Theta] * kappa1_times_lapse *
                (((double)GR_SPACEDIM + 1) +
                 m_params.kappa2 * ((double)GR_SPACEDIM - 1.0)) -
            Z_dot_d1lapse - state_cell_data[c_lapse] * m_cosmological_constant;

        rhs_cell_data[c_K] =
            advec_K +
            state_cell_data[c_lapse] *
                (ricci.scalar +
                 state_cell_data[c_K] *
                     (state_cell_data[c_K] - 2.0 * state_cell_data[c_Theta])) -
            kappa1_times_lapse * (double)GR_SPACEDIM * (1.0 + m_params.kappa2) *
                state_cell_data[c_Theta] -
            tr_covd2lapse -
            2.0 * state_cell_data[c_lapse] * (double)GR_SPACEDIM /
                ((double)GR_SPACEDIM - 1.0) * m_cosmological_constant;
    }

    // Gamma specific parts:
    //    Tensor<1, amrex::Real> Gammadot;
    auto d2_shift = m_deriv.diff2_sym_vector(ix, iy, iz, state, c_shift1);
    auto d1_K     = m_deriv.diff1_array_scalar(ix, iy, iz, state, c_K);
    auto d1_Theta = m_deriv.diff1_array_scalar(ix, iy, iz, state, c_Theta);

    FOR (i)
    {
        rhs_cell_data[c_Gamma1 + i] =
            (2.0 / (double)GR_SPACEDIM) *
                (divshift * (chris.contracted(i) +
                             2.0 * m_params.kappa3 * Z_over_chi[i]) -
                 2.0 * state_cell_data[c_lapse] * state_cell_data[c_K] *
                     Z_over_chi[i]) -
            2.0 * kappa1_times_lapse * Z_over_chi[i];

        FOR (j)
        {
            //            int idx1 = i + j + ((i * j != 0) ? 1 : 0);
            rhs_cell_data[c_Gamma1 + i] +=
                2.0 * h_UU(i, j) *
                    (state_cell_data[c_lapse] * d1_Theta(j) -
                     state_cell_data[c_Theta] * d1_lapse[j]) -
                2.0 * A_UU(i, j) * d1_lapse[j] -
                state_cell_data[c_lapse] *
                    ((2.0 * ((double)GR_SPACEDIM - 1.0) / (double)GR_SPACEDIM) *
                         h_UU(i, j) * d1_K(j) +
                     (double)GR_SPACEDIM * A_UU(i, j) * d1_chi[j] /
                         state_cell_data[c_chi]) -
                (chris.contracted(j) + 2.0 * m_params.kappa3 * Z_over_chi[j]) *
                    d1_shift[i][j];

            FOR (k)
            {
                //                int idx2 = j + k + ((j * k != 0) ? 1 : 0);
                rhs_cell_data[c_Gamma1 + i] +=
                    2.0 * state_cell_data[c_lapse] * chris.ULL(i, j, k) *

                        A_UU(j, k) +
                    h_UU(j, k) * d2_shift(i, SYMM_IDX(j, k)) +
                    (((double)GR_SPACEDIM - 2.0) / (double)GR_SPACEDIM) *
                        h_UU(i, j) * d2_shift(k, SYMM_IDX(j, k));
            }
        }
    }

    auto advec_Gamma =
        m_deriv.advec_vector(ix, iy, iz, state, shift_vector, c_Gamma1);
    FOR (i)
    {
        rhs_cell_data[c_Gamma1 + i] += advec_Gamma[i];
    }
}

template <class gauge_t, class deriv_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
CCZ4RHS<gauge_t, deriv_t>::calculate_A_ij_rhs_no_amrex_array(
    int ix, int iy, int iz, const amrex::Array4<amrex::Real> &rhs_state,
    const amrex::Array4<amrex::Real const> &state) const
{
    const amrex::CellData<amrex::Real> &rhs_cell_data =
        rhs_state.cellData(ix, iy, iz);
    const amrex::CellData<const amrex::Real> &state_cell_data =
        state.cellData(ix, iy, iz);

    const auto h_UU = CCZ4Geometry::compute_inverse_metric(state_cell_data);

    // hij derivatives
    Tensor<3, amrex::Real> d1_h =
        m_deriv.diff1_tensor(ix, iy, iz, state, c_h11);

    const auto chris = CCZ4Geometry::compute_christoffel(d1_h, h_UU);

    Tensor<1, amrex::Real> Z_over_chi;
    Tensor<1, amrex::Real> Z; // NOLINT(readability-identifier-length)

    if (m_formulation == USE_BSSN)
    {
        FOR (i)
            Z_over_chi[i] = 0.0;
    }
    else
    {

        FOR (i)
            Z_over_chi[i] =
                0.5 * (state_cell_data[c_Gamma1 + i] - chris.contracted[i]);
    }
    FOR (i)
        Z[i] = state_cell_data[c_chi] * Z_over_chi[i];

    auto ricci = CCZ4Geometry::compute_ricci_Z(
        ix, iy, iz, rhs_state, state, h_UU, chris, Z_over_chi, m_deriv);

    Tensor<2, amrex::Real> d1_shift =
        m_deriv.diff1_vector(ix, iy, iz, state, c_shift1);
    amrex::Real divshift = CCZ4Geometry::compute_divshift(d1_shift);

    auto d1_lapse = m_deriv.diff1_scalar(ix, iy, iz, state, c_lapse);
    auto d1_chi   = m_deriv.diff1_scalar(ix, iy, iz, state, c_chi);

    amrex::Real Z_dot_d1lapse = TensorAlgebra::compute_dot_product(Z, d1_lapse);
    amrex::Real dlapse_dot_dchi =
        TensorAlgebra::compute_dot_product(d1_lapse, d1_chi, h_UU);

    Tensor<2, amrex::Real> covdtilde2lapse;
    Tensor<2, amrex::Real> covd2lapse;
    auto d2_lapse = m_deriv.diff2_sym_scalar(ix, iy, iz, state, c_lapse);

    FOR (k, l)
    {
        int idx               = k + l + ((k * l != 0) ? 1 : 0);
        covdtilde2lapse[k][l] = d2_lapse(idx);
        FOR (m)
        {
            covdtilde2lapse[k][l] -= chris.ULL[m][k][l] * d1_lapse[m];
        }
        covd2lapse[k][l] =
            state_cell_data[c_chi] * covdtilde2lapse[k][l] +
            0.5 * (d1_lapse[k] * d1_chi[l] + d1_chi[k] * d1_lapse[l] -
                   state_cell_data[var_idx(c_h11, k, l)] * dlapse_dot_dchi);
    }

    Tensor<2, amrex::Real> A_UU =
        CCZ4Geometry::compute_A_UU(state_cell_data, h_UU);

    Tensor<2, amrex::Real> Adot_TF;
    FOR (i, j)
    {
        Adot_TF[i][j] = -covd2lapse[i][j] + state_cell_data[c_chi] *
                                                state_cell_data[c_lapse] *
                                                ricci.LL[i][j];
    }
    CCZ4Geometry::make_trace_free(Adot_TF, state_cell_data, h_UU);

    Tensor<1, amrex::Real> shift_vector;
    FOR (idir)
    {
        shift_vector[idir] = state(ix, iy, iz, c_shift1 + idir);
    }
    auto advec_A = m_deriv.advec_tensor(ix, iy, iz, state, shift_vector, c_A11);

    FOR (i, j)
    {
        rhs_cell_data[var_idx(c_A11, i, j)] =
            advec_A[i][j] + Adot_TF[i][j] +
            state_cell_data[var_idx(c_A11, i, j)] *
                (state_cell_data[c_lapse] *
                     (state_cell_data[c_K] - 2.0 * state_cell_data[c_Theta]) -
                 (2.0 / (double)GR_SPACEDIM) * divshift);
        FOR (k)
        {
            rhs_cell_data[var_idx(c_A11, i, j)] +=
                state_cell_data[var_idx(c_A11, k, i)] * d1_shift[k][j] +
                state_cell_data[var_idx(c_A11, k, j)] * d1_shift[k][i];
            FOR (l)
            {
                int idx = k + l + ((k * l != 0) ? 1 : 0);
                rhs_cell_data[var_idx(c_A11, i, j)] -=
                    2.0 * state_cell_data[c_lapse] * h_UU[k][l] *
                    state_cell_data[var_idx(c_A11, i, k)] *
                    state_cell_data[var_idx(c_A11, l, j)];
            }
        }
    }

    // Theta specific parts

    amrex::Real Aij_squared =
        CCZ4Geometry::compute_Aij_squared(state_cell_data, A_UU);

    amrex::Real advec_K =
        m_deriv.advection(ix, iy, iz, state, shift_vector, c_K);

    amrex::Real tr_covd2lapse = -((double)GR_SPACEDIM / 2.0) * dlapse_dot_dchi;
    FOR (i)
    {
        tr_covd2lapse -=
            state_cell_data[c_chi] * chris.contracted[i] * d1_lapse[i];
        FOR (j)
        {
            int idx = i + j + ((i * j != 0) ? 1 : 0);
            tr_covd2lapse +=
                h_UU[i][j] * (state_cell_data[c_chi] * d2_lapse(idx) +
                              d1_lapse[i] * d1_chi[j]);
        }
    }

    amrex::Real kappa1_times_lapse;
    if (m_params.covariantZ4)
    {
        kappa1_times_lapse = m_params.kappa1;
    }
    else
    {
        kappa1_times_lapse = m_params.kappa1 * state_cell_data[c_lapse];
    }

    if (m_formulation == USE_BSSN)
    {
        // ensure the Theta of CCZ4 remains at zero
        rhs_cell_data[c_Theta] = 0.0;
        // Use hamiltonian constraint to remove ricci.scalar for BSSN update
        rhs_cell_data[c_K] =
            advec_K +
            state_cell_data[c_lapse] *
                (Aij_squared + state_cell_data[c_K] * state_cell_data[c_K] /
                                   (double)GR_SPACEDIM) -
            tr_covd2lapse -
            2.0 * state_cell_data[c_lapse] * m_cosmological_constant /
                ((double)GR_SPACEDIM - 1.0);
    }
    else
    {
        amrex::Real advec_Theta =
            m_deriv.advection(ix, iy, iz, state, shift_vector, c_Theta);

        rhs_cell_data[c_Theta] =
            advec_Theta +
            0.5 * state_cell_data[c_lapse] *
                (ricci.scalar - Aij_squared +
                 (((double)GR_SPACEDIM - 1.0) / (double)GR_SPACEDIM) *
                     state_cell_data[c_K] * state_cell_data[c_K] -
                 2.0 * state_cell_data[c_Theta] * state_cell_data[c_K]) -
            0.5 * state_cell_data[c_Theta] * kappa1_times_lapse *
                (((double)GR_SPACEDIM + 1) +
                 m_params.kappa2 * ((double)GR_SPACEDIM - 1.0)) -
            Z_dot_d1lapse - state_cell_data[c_lapse] * m_cosmological_constant;

        rhs_cell_data[c_K] =
            advec_K +
            state_cell_data[c_lapse] *
                (ricci.scalar +
                 state_cell_data[c_K] *
                     (state_cell_data[c_K] - 2.0 * state_cell_data[c_Theta])) -
            kappa1_times_lapse * (double)GR_SPACEDIM * (1.0 + m_params.kappa2) *
                state_cell_data[c_Theta] -
            tr_covd2lapse -
            2.0 * state_cell_data[c_lapse] * (double)GR_SPACEDIM /
                ((double)GR_SPACEDIM - 1.0) * m_cosmological_constant;
    }

    // Gamma specific parts:
    Tensor<1, amrex::Real> Gammadot;
    auto d2_shift = m_deriv.diff2_sym_vector(ix, iy, iz, state, c_shift1);
    auto d1_K     = m_deriv.diff1_scalar(ix, iy, iz, state, c_K);
    auto d1_Theta = m_deriv.diff1_scalar(ix, iy, iz, state, c_Theta);

    FOR (i)
    {
        Gammadot[i] = (2.0 / (double)GR_SPACEDIM) *
                          (divshift * (chris.contracted[i] +
                                       2.0 * m_params.kappa3 * Z_over_chi[i]) -
                           2.0 * state_cell_data[c_lapse] *
                               state_cell_data[c_K] * Z_over_chi[i]) -
                      2.0 * kappa1_times_lapse * Z_over_chi[i];

        FOR (j)
        {
            int idx1 = i + j + ((i * j != 0) ? 1 : 0);
            Gammadot[i] +=
                2.0 * h_UU[i][j] *
                    (state_cell_data[c_lapse] * d1_Theta[j] -
                     state_cell_data[c_Theta] * d1_lapse[j]) -
                2.0 * A_UU[i][j] * d1_lapse[j] -
                state_cell_data[c_lapse] *
                    ((2.0 * ((double)GR_SPACEDIM - 1.0) / (double)GR_SPACEDIM) *
                         h_UU[i][j] * d1_K[j] +
                     (double)GR_SPACEDIM * A_UU[i][j] * d1_chi[j] /
                         state_cell_data[c_chi]) -
                (chris.contracted[j] + 2.0 * m_params.kappa3 * Z_over_chi[j]) *
                    d1_shift[i][j];

            FOR (k)
            {
                int idx2 = j + k + ((j * k != 0) ? 1 : 0);
                Gammadot[i] +=
                    2.0 * state_cell_data[c_lapse] * chris.ULL[i][j][k] *
                        A_UU[j][k] +
                    h_UU[j][k] * d2_shift(i, idx2) +
                    (((double)GR_SPACEDIM - 2.0) / (double)GR_SPACEDIM) *
                        h_UU[i][j] * d2_shift(k, idx2);
            }
        }
    }

    auto advec_Gamma =
        m_deriv.advec_vector(ix, iy, iz, state, shift_vector, c_Gamma1);
    FOR (i)
    {
        rhs_cell_data[c_Gamma1 + i] = advec_Gamma[i] + Gammadot[i];
    }
}

template <class gauge_t, class deriv_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
CCZ4RHS<gauge_t, deriv_t>::calculate_Theta_rhs(
    int ix, int iy, int iz, const amrex::Array4<amrex::Real> &rhs_state,
    const amrex::Array4<amrex::Real const> &state) const
{
    const amrex::CellData<amrex::Real> &rhs_cell_data =
        rhs_state.cellData(ix, iy, iz);
    const amrex::CellData<const amrex::Real> &state_cell_data =
        state.cellData(ix, iy, iz);

    auto d1_lapse = m_deriv.diff1_scalar(ix, iy, iz, state, c_lapse);
    auto d2_lapse = m_deriv.diff2_scalar(ix, iy, iz, state, c_lapse);

    auto d1_chi = m_deriv.diff1_scalar(ix, iy, iz, state, c_chi);

    auto d1_h = m_deriv.diff1_tensor(ix, iy, iz, state, c_h11);

    const auto h_UU = CCZ4Geometry::compute_inverse_metric(state_cell_data);

    amrex::Real dlapse_dot_dchi =
        TensorAlgebra::compute_dot_product(d1_lapse, d1_chi, h_UU);

    const auto chris          = CCZ4Geometry::compute_christoffel(d1_h, h_UU);
    amrex::Real tr_covd2lapse = -((double)GR_SPACEDIM / 2.0) * dlapse_dot_dchi;
    FOR (i)
    {
        tr_covd2lapse -=
            state_cell_data[c_chi] * chris.contracted(i) * d1_lapse[i];
        FOR (j)
        {
            tr_covd2lapse +=
                h_UU[i][j] * (state_cell_data[c_chi] * d2_lapse[i][j] +
                              d1_lapse[i] * d1_chi[j]);
        }
    }

    Tensor<1, amrex::Real> shift_vector;
    FOR (idir)
    {
        shift_vector[idir] = state(ix, iy, iz, c_shift1 + idir);
    }

    amrex::Real advec_K =
        m_deriv.advection(ix, iy, iz, state, shift_vector, c_K);

    amrex::Real Aij_squared =
        CCZ4Geometry::compute_Aij_squared(state_cell_data, h_UU);

    if (m_formulation == USE_BSSN)
    {
        // ensure the Theta of CCZ4 remains at zero
        rhs_cell_data[c_Theta] = 0.0;
        // Use hamiltonian constraint to remove ricci.scalar for BSSN update
        rhs_cell_data[c_K] =
            advec_K +
            state_cell_data[c_lapse] *
                (Aij_squared + state_cell_data[c_K] * state_cell_data[c_K] /
                                   (double)GR_SPACEDIM) -
            tr_covd2lapse -
            2.0 * state_cell_data[c_lapse] * m_cosmological_constant /
                ((double)GR_SPACEDIM - 1.0);
    }
    else
    {
        amrex::Real advec_Theta =
            m_deriv.advection(ix, iy, iz, state, shift_vector, c_Theta);

        // hij derivatives
        Tensor<1, amrex::Real> Z_over_chi;

        FOR (i)
            Z_over_chi[i] =
                0.5 * (state_cell_data[c_Gamma1 + i] - chris.contracted[i]);

        Tensor<1, amrex::Real> Z; // NOLINT(readability-identifier-length)
        FOR (i)
            Z[i] = state_cell_data[c_chi] * Z_over_chi[i];

        auto ricci = CCZ4Geometry::compute_ricci_Z(
            ix, iy, iz, rhs_state, state, h_UU, chris, Z_over_chi, m_deriv);

        amrex::Real kappa1_times_lapse;
        if (m_params.covariantZ4)
        {
            kappa1_times_lapse = m_params.kappa1;
        }
        else
        {
            kappa1_times_lapse = m_params.kappa1 * state_cell_data[c_lapse];
        }

        amrex::Real Z_dot_d1lapse =
            TensorAlgebra::compute_dot_product(Z, d1_lapse);

        rhs_cell_data[c_Theta] =
            advec_Theta +
            0.5 * state_cell_data[c_lapse] *
                (ricci.scalar - Aij_squared +
                 (((double)GR_SPACEDIM - 1.0) / (double)GR_SPACEDIM) *
                     state_cell_data[c_K] * state_cell_data[c_K] -
                 2.0 * state_cell_data[c_Theta] * state_cell_data[c_K]) -
            0.5 * state_cell_data[c_Theta] * kappa1_times_lapse *
                (((double)GR_SPACEDIM + 1) +
                 m_params.kappa2 * ((double)GR_SPACEDIM - 1.0)) -
            Z_dot_d1lapse - state_cell_data[c_lapse] * m_cosmological_constant;

        rhs_cell_data[c_K] =
            advec_K +
            state_cell_data[c_lapse] *
                (ricci.scalar +
                 state_cell_data[c_K] *
                     (state_cell_data[c_K] - 2.0 * state_cell_data[c_Theta])) -
            kappa1_times_lapse * (double)GR_SPACEDIM * (1.0 + m_params.kappa2) *
                state_cell_data[c_Theta] -
            tr_covd2lapse -
            2.0 * state_cell_data[c_lapse] * (double)GR_SPACEDIM /
                ((double)GR_SPACEDIM - 1.0) * m_cosmological_constant;
    }
}

template <class gauge_t, class deriv_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
CCZ4RHS<gauge_t, deriv_t>::calculate_Gamma_rhs(
    int ix, int iy, int iz, const amrex::Array4<amrex::Real> &rhs_state,
    const amrex::Array4<amrex::Real const> &state) const
{
    const amrex::CellData<amrex::Real> &rhs_cell_data =
        rhs_state.cellData(ix, iy, iz);
    const amrex::CellData<const amrex::Real> &state_cell_data =
        state.cellData(ix, iy, iz);

    Tensor<2, amrex::Real> d1_shift =
        m_deriv.diff1_vector(ix, iy, iz, state, c_shift1);
    amrex::Real divshift = CCZ4Geometry::compute_divshift(d1_shift);

    // hij derivatives
    Tensor<3, amrex::Real> d1_h =
        m_deriv.diff1_tensor(ix, iy, iz, state, c_h11);

    const auto h_UU = CCZ4Geometry::compute_inverse_metric(state_cell_data);
    Tensor<2, amrex::Real> A_UU =
        CCZ4Geometry::compute_A_UU(state_cell_data, h_UU);

    const auto chris = CCZ4Geometry::compute_christoffel(d1_h, h_UU);

    Tensor<1, amrex::Real> Z_over_chi;
    Tensor<1, amrex::Real> Z; // NOLINT(readability-identifier-length)

    FOR (i)
        Z_over_chi[i] =
            0.5 * (state_cell_data[c_Gamma1 + i] - chris.contracted[i]);

    Tensor<1, amrex::Real> Gammadot;

    auto d1_lapse = m_deriv.diff1_scalar(ix, iy, iz, state, c_lapse);
    auto d1_chi   = m_deriv.diff1_scalar(ix, iy, iz, state, c_chi);
    auto d1_K     = m_deriv.diff1_scalar(ix, iy, iz, state, c_K);
    auto d1_Theta = m_deriv.diff1_scalar(ix, iy, iz, state, c_Theta);
    auto d2_shift = m_deriv.diff2_vector(ix, iy, iz, state, c_shift1);

    amrex::Real kappa1_times_lapse;
    if (m_params.covariantZ4)
    {
        kappa1_times_lapse = m_params.kappa1;
    }
    else
    {
        kappa1_times_lapse = m_params.kappa1 * state_cell_data[c_lapse];
    }

    FOR (i)
    {
        Gammadot[i] = (2.0 / (double)GR_SPACEDIM) *
                          (divshift * (chris.contracted[i] +
                                       2.0 * m_params.kappa3 * Z_over_chi[i]) -
                           2.0 * state_cell_data[c_lapse] *
                               state_cell_data[c_K] * Z_over_chi[i]) -
                      2.0 * kappa1_times_lapse * Z_over_chi[i];

        FOR (j)
        {
            Gammadot[i] +=
                2.0 * h_UU[i][j] *
                    (state_cell_data[c_lapse] * d1_Theta[j] -
                     state_cell_data[c_Theta] * d1_lapse[j]) -
                2.0 * A_UU[i][j] * d1_lapse[j] -
                state_cell_data[c_lapse] *
                    ((2.0 * ((double)GR_SPACEDIM - 1.0) / (double)GR_SPACEDIM) *
                         h_UU[i][j] * d1_K[j] +
                     (double)GR_SPACEDIM * A_UU[i][j] * d1_chi[j] /
                         state_cell_data[c_chi]) -
                (chris.contracted[j] + 2.0 * m_params.kappa3 * Z_over_chi[j]) *
                    d1_shift[i][j];

            FOR (k)
            {
                Gammadot[i] +=
                    2.0 * state_cell_data[c_lapse] * chris.ULL[i][j][k] *
                        A_UU[j][k] +
                    h_UU[j][k] * d2_shift[i][j][k] +
                    (((double)GR_SPACEDIM - 2.0) / (double)GR_SPACEDIM) *
                        h_UU[i][j] * d2_shift[k][j][k];
            }
        }
    }
    Tensor<1, amrex::Real> shift_vector;
    FOR (idir)
    {
        shift_vector[idir] = state(ix, iy, iz, c_shift1 + idir);
    }

    auto advec_Gamma =
        m_deriv.advec_vector(ix, iy, iz, state, shift_vector, c_Gamma1);
    FOR (i)
    {
        rhs_cell_data[c_Gamma1 + i] = advec_Gamma[i] + Gammadot[i];
    }
}

template <class gauge_t, class deriv_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
CCZ4RHS<gauge_t, deriv_t>::apply_gauge_and_dissipation(
    int ix, int iy, int iz, const amrex::Array4<amrex::Real> &rhs_state,
    const amrex::Array4<amrex::Real const> &state) const
{
    const amrex::CellData<amrex::Real> &rhs_cell_data =
        rhs_state.cellData(ix, iy, iz);
    const amrex::CellData<const amrex::Real> &state_cell_data =
        state.cellData(ix, iy, iz);

    Tensor<1, amrex::Real> shift_vector;
    FOR (idir)
    {
        shift_vector[idir] = state(ix, iy, iz, c_shift1 + idir);
    }

    auto advec_lapse =
        m_deriv.advec_scalar(ix, iy, iz, state, shift_vector, c_lapse);

    auto advec_shift =
        m_deriv.advec_vector(ix, iy, iz, state, shift_vector, c_shift1);

    auto advec_B = m_deriv.advec_vector(ix, iy, iz, state, shift_vector, c_B1);

    auto advec_Gamma =
        m_deriv.advec_vector(ix, iy, iz, state, shift_vector, c_Gamma1);

    m_gauge.rhs_gauge(rhs_cell_data, state_cell_data, advec_lapse, advec_shift,
                      advec_B, advec_Gamma);

    m_deriv.add_dissipation(ix, iy, iz, rhs_cell_data, state, m_sigma,
                            NUM_CCZ4_VARS);
}
// NOLINTBEGIN(readability-function-cognitive-complexity)
template <class gauge_t, class deriv_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void

CCZ4RHS<gauge_t, deriv_t>::rhs_equation(const amrex::CellData<amrex::Real> &rhs,
                                        const CCZ4Vars &vars,
                                        const CCZ4D1Vars &d1,
                                        const CCZ4D2Vars &d2,
                                        const CCZ4AdvecVars &advec) const
{

    const auto h_UU  = CCZ4Geometry::compute_inverse_metric(vars);
    const auto chris = CCZ4Geometry::compute_christoffel(d1, h_UU);

    amrex::Array1D<amrex::Real, 0, 3> Z_over_chi;
    amrex::Array1D<amrex::Real, 0, 3>
        Z; // NOLINT(readability-identifier-length)

    if (m_formulation == USE_BSSN)
    {
        FOR (i)
            Z_over_chi(i) = 0.0;
    }
    else
    {

        FOR (i)
            Z_over_chi(i) = 0.5 * (vars.Gamma(i) - chris.contracted(i));
    }
    FOR (i)
      Z(i) = vars.chi() * Z_over_chi(i);

    auto ricci = CCZ4Geometry::compute_ricci_Z(vars, d1, d2.chi, d2.h, h_UU,
                                               chris, Z_over_chi);

    amrex::Real divshift = CCZ4Geometry::compute_divshift(d1);
    amrex::Real Z_dot_d1lapse =
        TensorAlgebra::compute_dot_product(Z, d1.lapse());
    amrex::Real dlapse_dot_dchi =
        TensorAlgebra::compute_dot_product(d1.lapse(), d1.chi(), h_UU);

    amrex::Array2D<amrex::Real, 0, 3, 0, 3> covdtilde2lapse;
    amrex::Array2D<amrex::Real, 0, 3, 0, 3> covd2lapse;
    FOR (k, l)
    {
        covdtilde2lapse(k, l) = d2.lapse[k][l];
        FOR (m)
        {
            covdtilde2lapse(k, l) -= chris.ULL(m, k, l) * d1.lapse(m);
        }
        covd2lapse(k, l) =
	  vars.chi() * covdtilde2lapse(k, l) +
            0.5 * (d1.lapse(k) * d1.chi(l) + d1.chi(k) * d1.lapse(l) -
                   vars.h(k, l) * dlapse_dot_dchi);
    }

    amrex::Real tr_covd2lapse = -((double)GR_SPACEDIM / 2.0) * dlapse_dot_dchi;
    FOR (i)
    {
      tr_covd2lapse -= vars.chi() * chris.contracted(i) * d1.lapse(i);
        FOR (j)
        {
	  tr_covd2lapse += h_UU(i, j) * (vars.chi() * d2.lapse[i][j] +
                                           d1.lapse(i) * d1.chi(j));
        }
    }

    amrex::Array2D<amrex::Real, 0, 3, 0, 3> A_UU = CCZ4Geometry::compute_A_UU(vars, h_UU);

    // A^{ij} A_{ij}
    amrex::Real Aij_squared =
        CCZ4Geometry::compute_Aij_squared_with_A_UU(vars, A_UU);
    rhs[c_chi] = advec.chi() + (2.0 / (double)GR_SPACEDIM) * vars.chi() *
                                   (vars.lapse() * vars.K() - divshift);

    FOR2_SYM(i, j)
    {

        rhs[VAR_IDX(c_h11, i, j)] =
            advec.h(i, j) - 2.0 * vars.lapse() * vars.A(i, j) -
            (2.0 / (double)GR_SPACEDIM) * vars.h(i, j) * divshift;
        FOR (k)
        {
            rhs[VAR_IDX(c_h11, i, j)] +=
                vars.h(k, i) * d1.shift(k, j) + vars.h(k, j) * d1.shift(k, i);
        }
    }

    amrex::Array2D<amrex::Real, 0, 3, 0, 3> Adot_TF;
    FOR (i, j)
    {
        Adot_TF(i, j) =
            -covd2lapse(i, j) + vars.chi() * vars.lapse() * ricci.LL(i, j);
    }
    CCZ4Geometry::make_trace_free(Adot_TF, vars, h_UU);

    FOR2_SYM(i, j)
    {

        rhs[VAR_IDX(c_A11, i, j)] =
            advec.A(i, j) + Adot_TF(i, j) +
            vars.A(i, j) * (vars.lapse() * (vars.K() - 2.0 * vars.Theta()) -
                            (2.0 / (double)GR_SPACEDIM) * divshift);
        FOR (k)
        {
            rhs[VAR_IDX(c_A11, i, j)] +=
                vars.A(k, i) * d1.shift(k, j) + vars.A(k, j) * d1.shift(k, i);
            FOR (l)
            {
                rhs[VAR_IDX(c_A11, i, j)] -= 2.0 * vars.lapse() * h_UU(k, l) *
                                             vars.A(i, k) * vars.A(l, j);
            }
        }
    }

    amrex::Real kappa1_times_lapse;
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
        // ensure the Theta of CCZ4 remains at zero
        rhs[c_Theta] = 0.0;
        // Use hamiltonian constraint to remove ricci.scalar for BSSN update
        rhs[c_K] = advec.K() +
                   vars.lapse() * (Aij_squared +
                                   vars.K() * vars.K() / (double)GR_SPACEDIM) -
                   tr_covd2lapse -
                   2.0 * vars.lapse() * m_cosmological_constant /
                       ((double)GR_SPACEDIM - 1.0);
    }
    else
    {
        rhs[c_Theta] =
            advec.Theta() +
            0.5 * vars.lapse() *
                (ricci.scalar - Aij_squared +
                 (((double)GR_SPACEDIM - 1.0) / (double)GR_SPACEDIM) *
                     vars.K() * vars.K() -
                 2.0 * vars.Theta() * vars.K()) -
            0.5 * vars.Theta() * kappa1_times_lapse *
                (((double)GR_SPACEDIM + 1) +
                 m_params.kappa2 * ((double)GR_SPACEDIM - 1.0)) -
            Z_dot_d1lapse - vars.lapse() * m_cosmological_constant;

        rhs[c_K] = advec.K() +
                   vars.lapse() * (ricci.scalar +
                                   vars.K() * (vars.K() - 2.0 * vars.Theta())) -
                   kappa1_times_lapse * (double)GR_SPACEDIM *
                       (1.0 + m_params.kappa2) * vars.Theta() -
                   tr_covd2lapse -
                   2.0 * vars.lapse() * (double)GR_SPACEDIM /
                       ((double)GR_SPACEDIM - 1.0) * m_cosmological_constant;
    }

    amrex::Array1D<amrex::Real, 0, 3> Gammadot;
    FOR (i)
    {

        Gammadot(i) = (2.0 / (double)GR_SPACEDIM) *
                          (divshift * (chris.contracted(i) +
                                       2.0 * m_params.kappa3 * Z_over_chi(i)) -
                           2.0 * vars.lapse() * vars.K() * Z_over_chi(i)) -
                      2.0 * kappa1_times_lapse * Z_over_chi(i);

        FOR (j)
        {
            Gammadot(i) +=
                2.0 * h_UU(i, j) *
                    (vars.lapse() * d1.Theta(j) - vars.Theta() * d1.lapse(j)) -
                2.0 * A_UU(i, j) * d1.lapse(j) -
                vars.lapse() *
                    ((2.0 * ((double)GR_SPACEDIM - 1.0) / (double)GR_SPACEDIM) *
                         h_UU(i, j) * d1.K(j) +
                     (double)GR_SPACEDIM * A_UU(i, j) * d1.chi(j) /
                         vars.chi()) -
                (chris.contracted(j) + 2.0 * m_params.kappa3 * Z_over_chi(j)) *
                    d1.shift(i, j);

            FOR (k)
            {
                Gammadot(i) +=
                    2.0 * vars.lapse() * chris.ULL(i, j, k) * A_UU(j, k) +
                    h_UU(j, k) * d2.shift[i][j][k] +
                    (((double)GR_SPACEDIM - 2.0) / (double)GR_SPACEDIM) *
                        h_UU(i, j) * d2.shift[k][j][k];
            }
        }
    }
    FOR (i)
    {
        rhs[c_Gamma1 + i] = advec.Gamma(i) + Gammadot(i);
    }

    m_gauge.rhs_gauge(rhs, vars, d1, d2, advec);
}

template <class gauge_t, class deriv_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
CCZ4RHS<gauge_t, deriv_t>::calculate_rhs(
    int ix, int iy, int iz, const amrex::Array4<amrex::Real> &rhs_state,
    const amrex::Array4<amrex::Real const> &state) const
{
    const amrex::CellData<amrex::Real> &rhs_cell_data =
        rhs_state.cellData(ix, iy, iz);
    const amrex::CellData<const amrex::Real> &state_cell_data =
        state.cellData(ix, iy, iz);

    Tensor<1, amrex::Real> shift_vector;
    FOR (idir)
    {
        shift_vector[idir] = state(ix, iy, iz, c_shift1 + idir);
    }

    // Tensor<2, amrex::Real> d1_shift =
    //     m_deriv.diff1_vector(ix, iy, iz, state, c_shift1);
    auto d1_shift = m_deriv.diff1_array_vector(ix, iy, iz, state, c_shift1);
    amrex::Real divshift = CCZ4Geometry::compute_divshift(d1_shift);

    amrex::Real advec_chi =
        m_deriv.advection(ix, iy, iz, state, shift_vector, c_chi);
    rhs_cell_data[c_chi] =
        advec_chi +
        (2.0 / (double)GR_SPACEDIM) * state_cell_data[c_chi] *
            (state_cell_data[c_lapse] * state_cell_data[c_K] - divshift);

    // Calculation of h_ij RHS

    const auto advec_h =
        m_deriv.advec_tensor(ix, iy, iz, state, shift_vector, c_h11);

    FOR (i, j)
    {
        int idx = SYMM_IDX(i, j);
        rhs_cell_data[c_h11 + idx] =
            advec_h[i][j] -
            2.0 * state_cell_data[c_lapse] * state_cell_data[c_A11 + idx] -
            (2.0 / (double)GR_SPACEDIM) * state_cell_data[c_h11 + idx] *
                divshift;

        FOR (k)
        {
            rhs_cell_data[c_h11 + idx] +=
                state_cell_data[var_idx(c_h11, k, i)] * d1_shift(k, j) +
                state_cell_data[var_idx(c_h11, k, j)] * d1_shift(k, i);
        }
    }
    /////////////////////////////////////////////////////////////////////////
    const amrex::Array2D<amrex::Real, 0, 3, 0, 3> h_UU =
        CCZ4Geometry::compute_inverse_metric_array(state_cell_data);

    // hij derivatives
    const auto d1_h = m_deriv.diff1_array_tensor(ix, iy, iz, state, c_h11);

    const chris_t_array chris = CCZ4Geometry::compute_christoffel(d1_h, h_UU);

    amrex::Array1D<amrex::Real, 0, 3> Z_over_chi;
    amrex::Array1D<amrex::Real, 0, 3>
        Z; // NOLINT(readability-identifier-length)

    if (m_formulation == USE_BSSN)
    {
        FOR (i)
            Z_over_chi(i) = 0.0;
    }
    else
    {

        FOR (i)
            Z_over_chi(i) =
                0.5 * (state_cell_data[c_Gamma1 + i] - chris.contracted(i));
    }
    FOR (i)
        Z(i) = state_cell_data[c_chi] * Z_over_chi(i);

    const auto ricci = CCZ4Geometry::compute_ricci_Z(
        ix, iy, iz, rhs_state, state, h_UU, chris, Z_over_chi, m_deriv);

    // auto d1_shift = m_deriv.diff1_array_vector(ix, iy, iz, state, c_shift1);
    // amrex::Real divshift = CCZ4Geometry::compute_divshift(d1_shift);

    const auto d1_lapse =
        m_deriv.diff1_array_scalar(ix, iy, iz, state, c_lapse);
    const auto d1_chi = m_deriv.diff1_array_scalar(ix, iy, iz, state, c_chi);

    const amrex::Real Z_dot_d1lapse =
        TensorAlgebra::compute_dot_product(Z, d1_lapse);
    const amrex::Real dlapse_dot_dchi =
        TensorAlgebra::compute_dot_product(d1_lapse, d1_chi, h_UU);

    amrex::Array2D<amrex::Real, 0, 3, 0, 3> covdtilde2lapse;
    amrex::Array2D<amrex::Real, 0, 3, 0, 3> covd2lapse;
    const auto d2_lapse =
        m_deriv.diff2_array_scalar(ix, iy, iz, state, c_lapse);

    FOR (k, l)
    {
        covdtilde2lapse(k, l) = d2_lapse(k, l);
        FOR (m)
        {
            covdtilde2lapse(k, l) -= chris.ULL(m, k, l) * d1_lapse(m);
        }
        covd2lapse(k, l) =
            state_cell_data[c_chi] * covdtilde2lapse(k, l) +
            0.5 * (d1_lapse(k) * d1_chi(l) + d1_chi(k) * d1_lapse(l) -
                   state_cell_data[var_idx(c_h11, k, l)] * dlapse_dot_dchi);
    }

    const auto A_UU = CCZ4Geometry::compute_A_UU(state_cell_data, h_UU);

    amrex::Array2D<amrex::Real, 0, 3, 0, 3> Adot_TF;
    FOR (i, j)
    {
        Adot_TF(i, j) = -covd2lapse(i, j) + state_cell_data[c_chi] *
                                                state_cell_data[c_lapse] *
                                                ricci.LL(i, j);
    }
    CCZ4Geometry::make_trace_free(Adot_TF, state_cell_data, h_UU);

    // Tensor<1, amrex::Real> shift_vector;
    // FOR (idir)
    // {
    //     shift_vector[idir] = state(ix, iy, iz, c_shift1 + idir);
    // }
    const auto advec_A =
        m_deriv.advec_tensor(ix, iy, iz, state, shift_vector, c_A11);

    FOR (i, j)
    {
        int idx = SYMM_IDX(i, j);
        rhs_cell_data[c_A11 + idx] =
            advec_A[i][j] + Adot_TF(i, j) +
            state_cell_data[c_A11 + idx] *
                (state_cell_data[c_lapse] *
                     (state_cell_data[c_K] - 2.0 * state_cell_data[c_Theta]) -
                 (2.0 / (double)GR_SPACEDIM) * divshift);
        FOR (k)
        {
            rhs_cell_data[c_A11 + idx] +=
                state_cell_data[var_idx(c_A11, k, i)] * d1_shift(k, j) +
                state_cell_data[var_idx(c_A11, k, j)] * d1_shift(k, i);
            FOR (l)
            {
                rhs_cell_data[c_A11 + idx] -=
                    2.0 * state_cell_data[c_lapse] * h_UU(k, l) *
                    state_cell_data[var_idx(c_A11, i, k)] *
                    state_cell_data[var_idx(c_A11, l, j)];
            }
        }
    }

    // Theta specific parts

    const amrex::Real Aij_squared =
        CCZ4Geometry::compute_Aij_squared(state_cell_data, A_UU);

    const amrex::Real advec_K =
        m_deriv.advection(ix, iy, iz, state, shift_vector, c_K);

    amrex::Real tr_covd2lapse = -((double)GR_SPACEDIM / 2.0) * dlapse_dot_dchi;
    FOR (i)
    {
        tr_covd2lapse -=
            state_cell_data[c_chi] * chris.contracted(i) * d1_lapse(i);
        FOR (j)
        {
            tr_covd2lapse +=
                h_UU(i, j) * (state_cell_data[c_chi] * d2_lapse(i, j) +
                              d1_lapse(i) * d1_chi(j));
        }
    }

    amrex::Real kappa1_times_lapse;
    if (m_params.covariantZ4)
    {
        kappa1_times_lapse = m_params.kappa1;
    }
    else
    {
        kappa1_times_lapse = m_params.kappa1 * state_cell_data[c_lapse];
    }

    if (m_formulation == USE_BSSN)
    {
        // ensure the Theta of CCZ4 remains at zero
        rhs_cell_data[c_Theta] = 0.0;
        // Use hamiltonian constraint to remove ricci.scalar for BSSN update
        rhs_cell_data[c_K] =
            advec_K +
            state_cell_data[c_lapse] *
                (Aij_squared + state_cell_data[c_K] * state_cell_data[c_K] /
                                   (double)GR_SPACEDIM) -
            tr_covd2lapse -
            2.0 * state_cell_data[c_lapse] * m_cosmological_constant /
                ((double)GR_SPACEDIM - 1.0);
    }
    else
    {
        amrex::Real advec_Theta =
            m_deriv.advection(ix, iy, iz, state, shift_vector, c_Theta);

        rhs_cell_data[c_Theta] =
            advec_Theta +
            0.5 * state_cell_data[c_lapse] *
                (ricci.scalar - Aij_squared +
                 (((double)GR_SPACEDIM - 1.0) / (double)GR_SPACEDIM) *
                     state_cell_data[c_K] * state_cell_data[c_K] -
                 2.0 * state_cell_data[c_Theta] * state_cell_data[c_K]) -
            0.5 * state_cell_data[c_Theta] * kappa1_times_lapse *
                (((double)GR_SPACEDIM + 1) +
                 m_params.kappa2 * ((double)GR_SPACEDIM - 1.0)) -
            Z_dot_d1lapse - state_cell_data[c_lapse] * m_cosmological_constant;

        rhs_cell_data[c_K] =
            advec_K +
            state_cell_data[c_lapse] *
                (ricci.scalar +
                 state_cell_data[c_K] *
                     (state_cell_data[c_K] - 2.0 * state_cell_data[c_Theta])) -
            kappa1_times_lapse * (double)GR_SPACEDIM * (1.0 + m_params.kappa2) *
                state_cell_data[c_Theta] -
            tr_covd2lapse -
            2.0 * state_cell_data[c_lapse] * (double)GR_SPACEDIM /
                ((double)GR_SPACEDIM - 1.0) * m_cosmological_constant;
    }

    // Gamma specific parts:
    //    Tensor<1, amrex::Real> Gammadot;
    const auto d2_shift =
        m_deriv.diff2_array_vector(ix, iy, iz, state, c_shift1);
    const auto d1_K = m_deriv.diff1_array_scalar(ix, iy, iz, state, c_K);
    const auto d1_Theta =
        m_deriv.diff1_array_scalar(ix, iy, iz, state, c_Theta);

    const auto advec_Gamma =
        m_deriv.advec_vector(ix, iy, iz, state, shift_vector, c_Gamma1);

    FOR (i)
    {
        rhs_cell_data[c_Gamma1 + i] =
            (2.0 / (double)GR_SPACEDIM) *
                (divshift * (chris.contracted(i) +
                             2.0 * m_params.kappa3 * Z_over_chi(i)) -
                 2.0 * state_cell_data[c_lapse] * state_cell_data[c_K] *
                     Z_over_chi(i)) -
            2.0 * kappa1_times_lapse * Z_over_chi(i) + advec_Gamma[i];

        FOR (j)
        {
            //            int idx1 = i + j + ((i * j != 0) ? 1 : 0);
            rhs_cell_data[c_Gamma1 + i] +=
                2.0 * h_UU(i, j) *
                    (state_cell_data[c_lapse] * d1_Theta(j) -
                     state_cell_data[c_Theta] * d1_lapse(j)) -
                2.0 * A_UU(i, j) * d1_lapse(j) -
                state_cell_data[c_lapse] *
                    ((2.0 * ((double)GR_SPACEDIM - 1.0) / (double)GR_SPACEDIM) *
                         h_UU(i, j) * d1_K(j) +
                     (double)GR_SPACEDIM * A_UU(i, j) * d1_chi(j) /
                         state_cell_data[c_chi]) -
                (chris.contracted(j) + 2.0 * m_params.kappa3 * Z_over_chi(j)) *
                    d1_shift(i, j);

            FOR (k)
            {
                //                int idx2 = j + k + ((j * k != 0) ? 1 : 0);
                rhs_cell_data[c_Gamma1 + i] +=
                    2.0 * state_cell_data[c_lapse] * chris.ULL(i, j, k) *
                        A_UU(j, k) +
                    h_UU(j, k) * d2_shift(i, j, k) +
                    (((double)GR_SPACEDIM - 2.0) / (double)GR_SPACEDIM) *
                        h_UU(i, j) * d2_shift(k, j, k);
            }
        }
    }

    //////////////////////////////////////////////
    const auto advec_lapse =
        m_deriv.advec_scalar(ix, iy, iz, state, shift_vector, c_lapse);

    const auto advec_shift =
        m_deriv.advec_vector(ix, iy, iz, state, shift_vector, c_shift1);

    const auto advec_B =
        m_deriv.advec_vector(ix, iy, iz, state, shift_vector, c_B1);

    // auto advec_Gamma =
    //     m_deriv.advec_vector(ix, iy, iz, state, shift_vector, c_Gamma1);

    m_gauge.rhs_gauge(rhs_cell_data, state_cell_data, advec_lapse, advec_shift,
                      advec_B, advec_Gamma);

    m_deriv.add_dissipation(ix, iy, iz, rhs_cell_data, state, m_sigma,
                            NUM_CCZ4_VARS);
}

// NOLINTEND(readability-function-cognitive-complexity)

#endif /* CCZ4RHS_IMPL_HPP_ */
