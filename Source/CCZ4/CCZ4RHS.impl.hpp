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
    int ix, int iy, int iz, const amrex::Array4<amrex::Real> &rhs,
    const amrex::Array4<amrex::Real const> &state) const
{
    const amrex::CellData<const amrex::Real> &state_cell_data =
        state.cellData(ix, iy, iz);
    ConstCCZ4Vars vars(state_cell_data);

    // Get the derivatives
    const CCZ4D1Vars d1(ix, iy, iz, state, m_deriv);
    const CCZ4D2Vars d2(ix, iy, iz, state, m_deriv);
    const CCZ4AdvecVars advec(ix, iy, iz, state, m_deriv);

    const amrex::CellData<const amrex::Real> &rhs_cell_data =
        rhs.cellData(ix, iy, iz);
    CCZ4Vars2 rhs_vars(rhs_cell_data);

    rhs_equation(rhs_vars, vars, d1, d2, advec);

    m_deriv.add_dissipation(ix, iy, iz, rhs_vars, state, m_sigma);
}

// NOLINTBEGIN(readability-function-cognitive-complexity)
template <class gauge_t, class deriv_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
CCZ4RHS<gauge_t, deriv_t>::rhs_equation(CCZ4Vars2 &rhs,
                                        const ConstCCZ4Vars &vars,
                                        const CCZ4D1Vars &d1,
                                        const CCZ4D2Vars &d2,
                                        const CCZ4AdvecVars &advec) const
{

    const auto h_UU  = CCZ4Geometry2::compute_inverse_metric(vars);
    const auto chris = TensorAlgebra::compute_christoffel(d1.h, h_UU);

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

    auto ricci = CCZ4Geometry2::compute_ricci_Z(vars, d1, d2.chi, d2.h, h_UU,
                                                chris, Z_over_chi);

    amrex::Real divshift      = TensorAlgebra::compute_trace(d1.shift);
    amrex::Real Z_dot_d1lapse = TensorAlgebra::compute_dot_product(Z, d1.lapse);
    amrex::Real dlapse_dot_dchi =
        TensorAlgebra::compute_dot_product(d1.lapse, d1.chi, h_UU);

    Tensor<2, amrex::Real> covdtilde2lapse;
    Tensor<2, amrex::Real> covd2lapse;
    FOR (k, l)
    {
        covdtilde2lapse[k][l] = d2.lapse[k][l];
        FOR (m)
        {
            covdtilde2lapse[k][l] -= chris.ULL[m][k][l] * d1.lapse[m];
        }
        covd2lapse[k][l] =
            vars.chi() * covdtilde2lapse[k][l] +
            0.5 * (d1.lapse[k] * d1.chi[l] + d1.chi[k] * d1.lapse[l] -
                   vars.h(k, l) * dlapse_dot_dchi);
    }

    amrex::Real tr_covd2lapse = -((double)GR_SPACEDIM / 2.0) * dlapse_dot_dchi;
    FOR (i)
    {
        tr_covd2lapse -= vars.chi() * chris.contracted[i] * d1.lapse[i];
        FOR (j)
        {
            tr_covd2lapse += h_UU[i][j] * (vars.chi() * d2.lapse[i][j] +
                                           d1.lapse[i] * d1.chi[j]);
        }
    }

    Tensor<2, amrex::Real> A_UU = CCZ4Geometry2::compute_A_UU(vars, h_UU);

    // A^{ij} A_{ij}
    amrex::Real Aij_squared = CCZ4Geometry2::compute_Aij_squared(vars, h_UU);
    amrex::Real rhs_chi = advec.chi + (2.0 / (double)GR_SPACEDIM) * vars.chi() *
                                          (vars.lapse() * vars.K() - divshift);
    rhs.store_chi(rhs_chi);

    Tensor<2, amrex::Real> rhs_h;
    FOR (i, j)
    {
        rhs_h[i][j] = advec.h[i][j] - 2.0 * vars.lapse() * vars.A(i, j) -
                      (2.0 / (double)GR_SPACEDIM) * vars.h(i, j) * divshift;
        FOR (k)
        {
            rhs_h[i][j] +=
                vars.h(k, i) * d1.shift[k][j] + vars.h(k, j) * d1.shift[k][i];
        }
    }
    rhs.store_h(rhs_h);

    Tensor<2, amrex::Real> Adot_TF;
    FOR (i, j)
    {
        Adot_TF[i][j] =
            -covd2lapse[i][j] + vars.chi() * vars.lapse() * ricci.LL[i][j];
    }
    CCZ4Geometry2::make_trace_free(Adot_TF, vars, h_UU);

    Tensor<2, amrex::Real> rhs_A;
    FOR (i, j)
    {
        rhs_A[i][j] =
            advec.A[i][j] + Adot_TF[i][j] +
            vars.A(i, j) * (vars.lapse() * (vars.K() - 2.0 * vars.Theta()) -
                            (2.0 / (double)GR_SPACEDIM) * divshift);
        FOR (k)
        {
            rhs_A[i][j] +=
                vars.A(k, i) * d1.shift[k][j] + vars.A(k, j) * d1.shift[k][i];
            FOR (l)
            {
                rhs_A[i][j] -= 2.0 * vars.lapse() * h_UU[k][l] * vars.A(i, k) *
                               vars.A(l, j);
            }
        }
    }
    rhs.store_A(rhs_A);

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
        rhs.store_Theta(0.0);
        // Use hamiltonian constraint to remove ricci.scalar for BSSN update
        amrex::Real rhs_K =
            advec.K +
            vars.lapse() *
                (Aij_squared + vars.K() * vars.K() / (double)GR_SPACEDIM) -
            tr_covd2lapse -
            2.0 * vars.lapse() * m_cosmological_constant /
                ((double)GR_SPACEDIM - 1.0);
        rhs.store_K(rhs_K);
    }
    else
    {
        amrex::Real rhs_Theta =
            advec.Theta +
            0.5 * vars.lapse() *
                (ricci.scalar - Aij_squared +
                 (((double)GR_SPACEDIM - 1.0) / (double)GR_SPACEDIM) *
                     vars.K() * vars.K() -
                 2.0 * vars.Theta() * vars.K()) -
            0.5 * vars.Theta() * kappa1_times_lapse *
                (((double)GR_SPACEDIM + 1) +
                 m_params.kappa2 * ((double)GR_SPACEDIM - 1.0)) -
            Z_dot_d1lapse - vars.lapse() * m_cosmological_constant;
        rhs.store_Theta(rhs_Theta);

        amrex::Real rhs_K =
            advec.K +
            vars.lapse() *
                (ricci.scalar + vars.K() * (vars.K() - 2.0 * vars.Theta())) -
            kappa1_times_lapse * (double)GR_SPACEDIM * (1.0 + m_params.kappa2) *
                vars.Theta() -
            tr_covd2lapse -
            2.0 * vars.lapse() * (double)GR_SPACEDIM /
                ((double)GR_SPACEDIM - 1.0) * m_cosmological_constant;
        rhs.store_K(rhs_K);
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
                    (vars.lapse() * d1.Theta[j] - vars.Theta() * d1.lapse[j]) -
                2.0 * A_UU[i][j] * d1.lapse[j] -
                vars.lapse() *
                    ((2.0 * ((double)GR_SPACEDIM - 1.0) / (double)GR_SPACEDIM) *
                         h_UU[i][j] * d1.K[j] +
                     (double)GR_SPACEDIM * A_UU[i][j] * d1.chi[j] /
                         vars.chi()) -
                (chris.contracted[j] + 2.0 * m_params.kappa3 * Z_over_chi[j]) *
                    d1.shift[i][j];

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

    Tensor<1, amrex::Real> rhs_Gamma;
    FOR (i)
    {
        rhs_Gamma[i] = advec.Gamma[i] + Gammadot[i];
    }
    rhs.store_Gamma(rhs_Gamma);

    m_gauge.rhs_gauge(rhs, vars, d1, d2, advec);
}
// NOLINTEND(readability-function-cognitive-complexity)

#endif /* CCZ4RHS_IMPL_HPP_ */
