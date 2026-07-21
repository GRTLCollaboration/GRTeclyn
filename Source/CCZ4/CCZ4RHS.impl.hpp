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
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
CCZ4RHS<gauge_t, deriv_t>::compute_chi_and_h_ij(
    int ix, int iy, int iz, const amrex::Array4<amrex::Real> &rhs_state,
    const amrex::Array4<amrex::Real const> &state) const
{
    const amrex::CellData<amrex::Real> &rhs_cell_data =
        rhs_state.cellData(ix, iy, iz);
    const amrex::CellData<const amrex::Real> &state_cell_data =
        state.cellData(ix, iy, iz);

    CCZ4Vars vars(state_cell_data);

    Tensor::Rank1 shift_vector({vars.shift(0), vars.shift(1), vars.shift(2)});

    auto d1_shift        = m_deriv.d1_vector(ix, iy, iz, state, c_shift1);
    amrex::Real divshift = CCZ4Geometry::compute_divshift(d1_shift);
    amrex::Real advec_chi =
        m_deriv.advection(ix, iy, iz, state, shift_vector, c_chi);
    rhs_cell_data[c_chi] = advec_chi + (2.0 / (double)GR_SPACEDIM) *
                                           vars.chi() *
                                           (vars.lapse() * vars.K() - divshift);

    // Calculation of h_ij RHS

    FOR2_SYM(i, j)
    {
        rhs_cell_data[sym_var_idx(c_h11, i, j)] =
            m_deriv.advection(ix, iy, iz, state, shift_vector,
                              sym_var_idx(c_h11, i, j)) -
            2.0 * vars.lapse() * vars.A(i, j) -
            (2.0 / (double)GR_SPACEDIM) * vars.h(i, j) * divshift;

        FOR (k)
        {
            rhs_cell_data[sym_var_idx(c_h11, i, j)] +=
                vars.h(k, i) * d1_shift(k, j) + vars.h(k, j) * d1_shift(k, i);
        }
    }
}

template <class gauge_t, class deriv_t>
template <int formulation, int use_covariant_Z4>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
CCZ4RHS<gauge_t, deriv_t>::compute_A_ij_and_Theta_and_Gamma(
    int ix, int iy, int iz, const amrex::Array4<amrex::Real> &rhs_state,
    const amrex::Array4<amrex::Real const> &state) const
{
    const amrex::CellData<amrex::Real> &rhs_cell_data =
        rhs_state.cellData(ix, iy, iz);
    const amrex::CellData<const amrex::Real> &state_cell_data =
        state.cellData(ix, iy, iz);

    CCZ4Vars vars(state_cell_data);

    const auto h_UU = CCZ4Geometry::compute_inverse_metric_test(vars);

    // hij derivatives
    auto d1_h        = m_deriv.d1_sym_tensor(ix, iy, iz, state, c_h11);
    const auto chris = CCZ4Geometry::compute_christoffel(d1_h, h_UU);

    Tensor::Rank1 Z_over_chi;
    Tensor::Rank1 Z; // NOLINT(readability-identifier-length)

    if constexpr (formulation == USE_BSSN)
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

    // Gamma derivatives
    auto d1_Gamma = m_deriv.d1_vector(ix, iy, iz, state, c_Gamma1);

    // hij derivatives
    auto d2_h = m_deriv.d2_tensor(ix, iy, iz, state, c_h11);

    // chi derivatives
    auto d1_chi = m_deriv.d1_scalar(ix, iy, iz, state, c_chi);
    auto d2_chi = m_deriv.d2_scalar(ix, iy, iz, state, c_chi);

    auto ricci = CCZ4Geometry::compute_ricci_Z(
        vars, d1_chi, d1_Gamma, d1_h, d2_h, d2_chi, h_UU, chris, Z_over_chi);

    auto d1_shift        = m_deriv.d1_vector(ix, iy, iz, state, c_shift1);
    amrex::Real divshift = CCZ4Geometry::compute_divshift(d1_shift);

    auto d1_lapse = m_deriv.d1_scalar(ix, iy, iz, state, c_lapse);

    amrex::Real Z_dot_d1lapse = TensorAlgebra::compute_dot_product(Z, d1_lapse);
    amrex::Real dlapse_dot_dchi =
        TensorAlgebra::compute_dot_product(d1_lapse, d1_chi, h_UU);

    Tensor::Rank2 covdtilde2lapse{};
    Tensor::Rank2 covd2lapse{};
    auto d2_lapse = m_deriv.d2_scalar(ix, iy, iz, state, c_lapse);

    FOR (k, l)
    {
        covdtilde2lapse(k, l) = d2_lapse(k, l);
        FOR (m)
        {
            covdtilde2lapse(k, l) -= chris.ULL(m, k, l) * d1_lapse(m);
        }
        covd2lapse(k, l) =
            vars.chi() * covdtilde2lapse(k, l) +
            0.5 * (d1_lapse(k) * d1_chi(l) + d1_chi(k) * d1_lapse(l) -
                   vars.h(k, l) * dlapse_dot_dchi);
    }

    auto A_UU = CCZ4Geometry::compute_A_UU(state_cell_data, h_UU);

    Tensor::Rank2 Adot_TF;
    FOR (i, j)
    {
        Adot_TF(i, j) =
            -covd2lapse(i, j) + vars.chi() * vars.lapse() * ricci.LL(i, j);
    }
    CCZ4Geometry::make_trace_free(Adot_TF, state_cell_data, h_UU);

    Tensor::Rank1 shift_vector({vars.shift(0), vars.shift(1), vars.shift(2)});

    FOR2_SYM(i, j)
    {
        rhs_cell_data[sym_var_idx(c_A11, i, j)] =
            m_deriv.advection(ix, iy, iz, state, shift_vector,
                              sym_var_idx(c_A11, i, j)) +
            Adot_TF(i, j) +
            vars.A(i, j) * (vars.lapse() * (vars.K() - 2.0 * vars.Theta()) -
                            (2.0 / (double)GR_SPACEDIM) * divshift);
        FOR (k)
        {
            rhs_cell_data[sym_var_idx(c_A11, i, j)] +=
                vars.A(k, i) * d1_shift(k, j) + vars.A(k, j) * d1_shift(k, i);
            FOR (l)
            {
                rhs_cell_data[sym_var_idx(c_A11, i, j)] -=
                    2.0 * vars.lapse() * h_UU(k, l) * vars.A(i, k) *
                    vars.A(l, j);
            }
        }
    }

    // Theta specific parts

    amrex::Real Aij_squared =
        CCZ4Geometry::compute_Aij_squared_with_A_UU(state_cell_data, A_UU);

    amrex::Real advec_K =
        m_deriv.advection(ix, iy, iz, state, shift_vector, c_K);

    amrex::Real tr_covd2lapse = -((double)GR_SPACEDIM / 2.0) * dlapse_dot_dchi;
    FOR (i)
    {
        tr_covd2lapse -= vars.chi() * chris.contracted(i) * d1_lapse(i);
        FOR (j)
        {
            tr_covd2lapse += h_UU(i, j) * (vars.chi() * d2_lapse(i, j) +
                                           d1_lapse(i) * d1_chi(j));
        }
    }

    amrex::Real kappa1_times_lapse;
    if constexpr (use_covariant_Z4)
    {
        kappa1_times_lapse = m_params.kappa1;
    }
    else
    {
        kappa1_times_lapse = m_params.kappa1 * vars.lapse();
    }

    if constexpr (formulation == USE_BSSN)
    {
        // ensure the Theta of CCZ4 remains at zero
        rhs_cell_data[c_Theta] = 0.0;
        // Use hamiltonian constraint to remove ricci.scalar for BSSN update
        rhs_cell_data[c_K] =
            advec_K +
            vars.lapse() *
                (Aij_squared + vars.K() * vars.K() / (double)GR_SPACEDIM) -
            tr_covd2lapse -
            2.0 * vars.lapse() * m_cosmological_constant /
                ((double)GR_SPACEDIM - 1.0);
    }
    else
    {
        amrex::Real advec_Theta =
            m_deriv.advection(ix, iy, iz, state, shift_vector, c_Theta);

        rhs_cell_data[c_Theta] =
            advec_Theta +
            0.5 * vars.lapse() *
                (ricci.scalar - Aij_squared +
                 (((double)GR_SPACEDIM - 1.0) / (double)GR_SPACEDIM) *
                     vars.K() * vars.K() -
                 2.0 * vars.Theta() * vars.K()) -
            0.5 * vars.Theta() * kappa1_times_lapse *
                (((double)GR_SPACEDIM + 1) +
                 m_params.kappa2 * ((double)GR_SPACEDIM - 1.0)) -
            Z_dot_d1lapse - vars.lapse() * m_cosmological_constant;

        rhs_cell_data[c_K] =
            advec_K +
            vars.lapse() *
                (ricci.scalar + vars.K() * (vars.K() - 2.0 * vars.Theta())) -
            kappa1_times_lapse * (double)GR_SPACEDIM * (1.0 + m_params.kappa2) *
                vars.Theta() -
            tr_covd2lapse -
            2.0 * vars.lapse() * (double)GR_SPACEDIM /
                ((double)GR_SPACEDIM - 1.0) * m_cosmological_constant;
    }

    // Gamma specific parts:
    auto d2_shift = m_deriv.d2_vector(ix, iy, iz, state, c_shift1);
    auto d1_K     = m_deriv.d1_scalar(ix, iy, iz, state, c_K);
    auto d1_Theta = m_deriv.d1_scalar(ix, iy, iz, state, c_Theta);

    FOR (i)
    {
        rhs_cell_data[c_Gamma1 + i] =
            (2.0 / (double)GR_SPACEDIM) *
                (divshift * (chris.contracted(i) +
                             2.0 * m_params.kappa3 * Z_over_chi(i)) -
                 2.0 * vars.lapse() * vars.K() * Z_over_chi(i)) -
            2.0 * kappa1_times_lapse * Z_over_chi(i);

        FOR (j)
        {
            //            int idx1 = i + j + ((i * j != 0) ? 1 : 0);
            rhs_cell_data[c_Gamma1 + i] +=
                2.0 * h_UU(i, j) *
                    (vars.lapse() * d1_Theta(j) - vars.Theta() * d1_lapse(j)) -
                2.0 * A_UU(i, j) * d1_lapse(j) -
                vars.lapse() *
                    ((2.0 * ((double)GR_SPACEDIM - 1.0) / (double)GR_SPACEDIM) *
                         h_UU(i, j) * d1_K(j) +
                     (double)GR_SPACEDIM * A_UU(i, j) * d1_chi(j) /
                         vars.chi()) -
                (chris.contracted(j) + 2.0 * m_params.kappa3 * Z_over_chi(j)) *
                    d1_shift(i, j);

            FOR (k)
            {
                //                int idx2 = j + k + ((j * k != 0) ? 1 : 0);
                rhs_cell_data[c_Gamma1 + i] +=
                    2.0 * vars.lapse() * chris.ULL(i, j, k) *

                        A_UU(j, k) +
                    h_UU(j, k) * d2_shift(i, j, k) +
                    (((double)GR_SPACEDIM - 2.0) / (double)GR_SPACEDIM) *
                        h_UU(i, j) * d2_shift(k, j, k);
            }
        }
    }

    auto advec_Gamma =
        m_deriv.advec_vector(ix, iy, iz, state, shift_vector, c_Gamma1);
    FOR (i)
    {
        rhs_cell_data[c_Gamma1 + i] += advec_Gamma(i);
    }
}

template <class gauge_t, class deriv_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void CCZ4RHS<gauge_t, deriv_t>::apply_gauge(
    int ix, int iy, int iz, const amrex::Array4<amrex::Real> &rhs_state,
    const amrex::Array4<amrex::Real const> &state) const
{
    const amrex::CellData<amrex::Real> &rhs_cell_data =
        rhs_state.cellData(ix, iy, iz);
    const amrex::CellData<const amrex::Real> &state_cell_data =
        state.cellData(ix, iy, iz);

    CCZ4Vars vars(state_cell_data);

    Tensor::Rank1 shift_vector({vars.shift(0), vars.shift(1), vars.shift(2)});

    auto advec_lapse =
        m_deriv.advec_scalar(ix, iy, iz, state, shift_vector, c_lapse);

    auto advec_shift =
        m_deriv.advec_vector(ix, iy, iz, state, shift_vector, c_shift1);

    auto advec_B = m_deriv.advec_vector(ix, iy, iz, state, shift_vector, c_B1);

    auto advec_Gamma =
        m_deriv.advec_vector(ix, iy, iz, state, shift_vector, c_Gamma1);

    m_gauge.rhs_gauge(rhs_cell_data, vars, advec_lapse, advec_shift, advec_B,
                      advec_Gamma);
}

template <class gauge_t, class deriv_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
CCZ4RHS<gauge_t, deriv_t>::apply_dissipation(
    int ix, int iy, int iz, const amrex::Array4<amrex::Real> &rhs_state,
    const amrex::Array4<amrex::Real const> &state) const
{
    const amrex::CellData<amrex::Real> &rhs_cell_data =
        rhs_state.cellData(ix, iy, iz);

    m_deriv.add_dissipation(ix, iy, iz, rhs_cell_data, state, m_sigma,
                            NUM_CCZ4_VARS);
}
// NOLINTEND(readability-function-cognitive-complexity)

#endif /* CCZ4RHS_IMPL_HPP_ */
