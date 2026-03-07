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

    // Get the derivatives
    const CCZ4D1Vars d1(ix, iy, iz, state, m_deriv);
    const CCZ4D2Vars d2(ix, iy, iz, state, m_deriv);
    const CCZ4AdvecVars advec(ix, iy, iz, state, m_deriv);

    const amrex::CellData<amrex::Real> &rhs_cell_data =
        rhs_state.cellData(ix, iy, iz);
    // CCZ4Vars rhs(rhs_cell_data);

    rhs_equation(rhs_cell_data, vars, d1, d2, advec);

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
            Z_over_chi[i] = 0.5 * (vars.Gamma(i) - chris.contracted[i]);
    }
    FOR (i)
        Z[i] = vars.chi() * Z_over_chi[i];

    auto ricci = CCZ4Geometry::compute_ricci_Z(vars, d1, d2.chi, d2.h, h_UU,
                                               chris, Z_over_chi);

    amrex::Real divshift = CCZ4Geometry::compute_divshift(d1);
    amrex::Real Z_dot_d1lapse =
        TensorAlgebra::compute_dot_product(Z, d1.lapse());
    amrex::Real dlapse_dot_dchi =
        TensorAlgebra::compute_dot_product(d1.lapse(), d1.chi(), h_UU);

    Tensor<2, amrex::Real> covdtilde2lapse;
    Tensor<2, amrex::Real> covd2lapse;
    FOR (k, l)
    {
        covdtilde2lapse[k][l] = d2.lapse[k][l];
        FOR (m)
        {
            covdtilde2lapse[k][l] -= chris.ULL[m][k][l] * d1.lapse(m);
        }
        covd2lapse[k][l] =
            vars.chi() * covdtilde2lapse[k][l] +
            0.5 * (d1.lapse(k) * d1.chi(l) + d1.chi(k) * d1.lapse(l) -
                   vars.h(k, l) * dlapse_dot_dchi);
    }

    amrex::Real tr_covd2lapse = -((double)GR_SPACEDIM / 2.0) * dlapse_dot_dchi;
    FOR (i)
    {
        tr_covd2lapse -= vars.chi() * chris.contracted[i] * d1.lapse(i);
        FOR (j)
        {
            tr_covd2lapse += h_UU[i][j] * (vars.chi() * d2.lapse[i][j] +
                                           d1.lapse(i) * d1.chi(j));
        }
    }

    Tensor<2, amrex::Real> A_UU = CCZ4Geometry::compute_A_UU(vars, h_UU);

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

    Tensor<2, amrex::Real> Adot_TF;
    FOR (i, j)
    {
        Adot_TF[i][j] =
            -covd2lapse[i][j] + vars.chi() * vars.lapse() * ricci.LL[i][j];
    }
    CCZ4Geometry::make_trace_free(Adot_TF, vars, h_UU);

    FOR2_SYM(i, j)
    {
        rhs[VAR_IDX(c_A11, i, j)] =
            advec.A(i, j) + Adot_TF[i][j] +
            vars.A(i, j) * (vars.lapse() * (vars.K() - 2.0 * vars.Theta()) -
                            (2.0 / (double)GR_SPACEDIM) * divshift);
        FOR (k)
        {
            rhs[VAR_IDX(c_A11, i, j)] +=
                vars.A(k, i) * d1.shift(k, j) + vars.A(k, j) * d1.shift(k, i);
            FOR (l)
            {
                rhs[VAR_IDX(c_A11, i, j)] -= 2.0 * vars.lapse() * h_UU[k][l] *
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

    Tensor<1, amrex::Real> Gammadot;
    FOR (i)
    {
        Gammadot[i] = (2.0 / (double)GR_SPACEDIM) *
                          (divshift * (chris.contracted[i] +
                                       2.0 * m_params.kappa3 * Z_over_chi[i]) -
                           2.0 * vars.lapse() * vars.K() * Z_over_chi[i]) -
                      2.0 * kappa1_times_lapse * Z_over_chi[i];

        FOR (j)
        {
            Gammadot[i] +=
                2.0 * h_UU[i][j] *
                    (vars.lapse() * d1.Theta(j) - vars.Theta() * d1.lapse(j)) -
                2.0 * A_UU[i][j] * d1.lapse(j) -
                vars.lapse() *
                    ((2.0 * ((double)GR_SPACEDIM - 1.0) / (double)GR_SPACEDIM) *
                         h_UU[i][j] * d1.K(j) +
                     (double)GR_SPACEDIM * A_UU[i][j] * d1.chi(j) /
                         vars.chi()) -
                (chris.contracted[j] + 2.0 * m_params.kappa3 * Z_over_chi[j]) *
                    d1.shift(i, j);

            FOR (k)
            {
                Gammadot[i] +=
                    2.0 * vars.lapse() * chris.ULL[i][j][k] * A_UU[j][k] +
                    h_UU[j][k] * d2.shift[i][j][k] +
                    (((double)GR_SPACEDIM - 2.0) / (double)GR_SPACEDIM) *
                        h_UU[i][j] * d2.shift[k][j][k];
            }
        }
    }

    FOR (i)
    {
        rhs[c_Gamma1 + i] = advec.Gamma(i) + Gammadot[i];
    }

    m_gauge.rhs_gauge(rhs, vars, d1, d2, advec);
}
// NOLINTEND(readability-function-cognitive-complexity)

#endif /* CCZ4RHS_IMPL_HPP_ */
