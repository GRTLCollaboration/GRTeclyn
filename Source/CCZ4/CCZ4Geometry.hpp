/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// This file calculates CCZ4 geometric quantities (or a similar 3+1 split).
#ifndef CCZ4GEOMETRY_HPP_
#define CCZ4GEOMETRY_HPP_

#include "CCZ4D1Vars.hpp"
#include "CCZ4Vars.hpp"
#include "DimensionDefinitions.hpp"
#include "TensorAlgebra.hpp"

//! A structure for the decomposed elements of the Energy Momentum Tensor in
//! 3+1D
struct emtensor_t
{
    Tensor<2, amrex::Real> S; //!< S_ij = T_ij
    Tensor<1, amrex::Real> j; //!< j_i = T_ia_n^a
    amrex::Real trS;          //!< trS = S^i_i
    amrex::Real rho;          //!< rho = T_ab n^a n^b
};

struct ricci_t
{
    Tensor<2, amrex::Real> LL; // Ricci with two indices down
    amrex::Real scalar{};      // Ricci scalar
};

namespace CCZ4Geometry
{
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_metric_determinant(const CCZ4Vars &vars)
{
    amrex::Real det = vars.h(0, 0) * vars.h(1, 1) * vars.h(2, 2) +
                      2 * vars.h(0, 1) * vars.h(0, 2) * vars.h(1, 2) -
                      vars.h(0, 0) * vars.h(1, 2) * vars.h(1, 2) -
                      vars.h(1, 1) * vars.h(0, 2) * vars.h(0, 2) -
                      vars.h(2, 2) * vars.h(0, 1) * vars.h(0, 1);

    return det;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor<2, amrex::Real>
compute_inverse_metric(const CCZ4Vars &vars)
{
    amrex::Real det_h         = compute_metric_determinant(vars);
    amrex::Real det_h_inverse = 1. / det_h;
    Tensor<2, amrex::Real> h_UU;
    h_UU[0][0] = (vars.h(1, 1) * vars.h(2, 2) - vars.h(1, 2) * vars.h(2, 1)) *
                 det_h_inverse;
    h_UU[0][1] = (vars.h(2, 0) * vars.h(1, 2) - vars.h(1, 0) * vars.h(2, 2)) *
                 det_h_inverse;
    h_UU[0][2] = (vars.h(1, 0) * vars.h(2, 1) - vars.h(2, 0) * vars.h(1, 1)) *
                 det_h_inverse;
    h_UU[1][1] = (vars.h(0, 0) * vars.h(2, 2) - vars.h(2, 0) * vars.h(0, 2)) *
                 det_h_inverse;
    h_UU[1][2] = (vars.h(0, 1) * vars.h(2, 0) - vars.h(0, 0) * vars.h(2, 1)) *
                 det_h_inverse;
    h_UU[2][2] = (vars.h(0, 0) * vars.h(1, 1) - vars.h(0, 1) * vars.h(1, 0)) *
                 det_h_inverse;
    h_UU[1][0] = h_UU[0][1];
    h_UU[2][0] = h_UU[0][2];
    h_UU[2][1] = h_UU[1][2];

    return h_UU;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_trace_A(const CCZ4Vars &vars)
{
    Tensor<2, amrex::Real> inverse_metric = compute_inverse_metric(vars);
    amrex::Real trace_A                   = 0.0;
    FOR (i, j)
    {
        trace_A += inverse_metric[i][j] * vars.A(i, j);
    }
    return trace_A;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_trace_A(const CCZ4Vars &vars,
                const Tensor<2, amrex::Real> &inverse_metric)
{
    amrex::Real trace_A = 0.0;
    FOR (i, j)
    {
        trace_A += inverse_metric[i][j] * vars.A(i, j);
    }
    return trace_A;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor<2, amrex::Real>
compute_A_UU(const CCZ4Vars &vars, const Tensor<2, amrex::Real> &inverse_metric)
{
    Tensor<2, amrex::Real> A_UU;
    FOR (i, j)
    {
        A_UU[i][j] = 0.0;
        FOR (k, l)
        {
            A_UU[i][j] +=
                inverse_metric[i][k] * inverse_metric[j][l] * vars.A(k, l);
        }
    }
    return A_UU;
}

// This is A_ij A^ij
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_Aij_squared(const CCZ4Vars &vars)
{
    Tensor<2, amrex::Real> inverse_metric = compute_inverse_metric(vars);
    amrex::Real Aij_squared               = 0.0;
    FOR (i, j, k, l)
    {
        Aij_squared += inverse_metric[i][k] * inverse_metric[j][l] *
                       vars.A(i, j) * vars.A(k, l);
    }
    return Aij_squared;
}

// This is A_ij A^ij
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_Aij_squared(const CCZ4Vars &vars,
                    const Tensor<2, amrex::Real> &inverse_metric)
{
    amrex::Real Aij_squared = 0.0;
    FOR (i, j, k, l)
    {
        Aij_squared += inverse_metric[i][k] * inverse_metric[j][l] *
                       vars.A(i, j) * vars.A(k, l);
    }
    return Aij_squared;
}

/// Computes the conformal christoffel symbol
AMREX_GPU_DEVICE AMREX_FORCE_INLINE chris_t
compute_christoffel(const CCZ4D1Vars &d1, const Tensor<2, amrex::Real> &h_UU)
{
    chris_t out{};

    FOR (i, j, k)
    {
        out.LLL[i][j][k] =
            0.5 * (d1.h(j, i, k) + d1.h(k, i, j) - d1.h(j, k, i));
    }
    FOR (i, j, k)
    {
        out.ULL[i][j][k] = 0;
        FOR (l)
        {
            out.ULL[i][j][k] += h_UU[i][l] * out.LLL[l][j][k];
        }
    }
    FOR (i)
    {
        out.contracted[i] = 0;
        FOR (j, k)
        {
            out.contracted[i] += h_UU[j][k] * out.ULL[i][j][k];
        }
    }

    return out;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor<3, amrex::Real>
compute_phys_chris(const Tensor<1, amrex::Real> &d1_chi, const CCZ4Vars &vars,
                   const Tensor<2, amrex::Real> &h_UU,
                   const Tensor<3, amrex::Real> &chris_ULL)
{
    using namespace TensorAlgebra;
    Tensor<3, amrex::Real> chris_phys;
    FOR (i, j, k)
    {
        chris_phys[i][j][k] =
            chris_ULL[i][j][k] -
            0.5 / vars.chi() *
                (delta(i, k) * d1_chi[j] + delta(i, j) * d1_chi[k]);
        FOR (m)
        {
            chris_phys[i][j][k] +=
                0.5 / vars.chi() * vars.h(j, k) * h_UU[i][m] * d1_chi[m];
        }
    }
    return chris_phys;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_divshift(const CCZ4D1Vars &d1)
{
    amrex::Real divshift = 0.;
    FOR (i)
        divshift += d1.shift(i, i);
    return divshift;
}

/// Removes the trace of a 2-Tensor with lower indices given a metric and an
/// inverse metric.

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
make_trace_free(Tensor<2, amrex::Real> &tensor_LL, const CCZ4Vars vars,
                const Tensor<2, amrex::Real> &inverse_metric)
{
    auto trace = TensorAlgebra::compute_trace(tensor_LL, inverse_metric);
    double one_over_gr_spacedim = 1. / ((double)GR_SPACEDIM);
    FOR (i, j)
    {
        tensor_LL[i][j] -= one_over_gr_spacedim * vars.h(i, j) * trace;
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_z_terms(const int i, const int j,
                const Tensor<1, amrex::Real> &Z_over_chi, const CCZ4Vars &vars,
                const Tensor<1, amrex::Real> &d1_chi)
{
    amrex::Real out = 0.;
    FOR (k)
    {
        out += Z_over_chi[k] *
               (vars.h(i, k) * d1_chi[j] + vars.h(j, k) * d1_chi[i] -
                vars.h(i, j) * d1_chi[k]);
    }
    return out;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE ricci_t compute_ricci_Z(
    const CCZ4Vars &vars, const CCZ4D1Vars &d1,
    const Tensor<2, amrex::Real> &d2_chi, const Tensor<4, amrex::Real> &d2_h,
    const Tensor<2, amrex::Real> &h_UU, const chris_t &chris,
    const Tensor<1, amrex::Real> &Z_over_chi)
{
    ricci_t out;

    Tensor<2, amrex::Real> covdtilde2chi;
    FOR (k, l)
    {
        covdtilde2chi[k][l] = d2_chi[k][l];
        FOR (m)
        {
            covdtilde2chi[k][l] -= chris.ULL[m][k][l] * d1.chi(m);
        }
    }

    Tensor<3, amrex::Real> chris_LLU = {0.};
    amrex::Real boxtildechi          = 0.;
    amrex::Real dchi_dot_dchi        = 0.;
    FOR (i, j)
    {
        boxtildechi   += covdtilde2chi[i][j] * h_UU[i][j];
        dchi_dot_dchi += d1.chi(i) * d1.chi(j) * h_UU[i][j];
        FOR (k, l)
        {
            chris_LLU[i][j][k] += h_UU[k][l] * chris.LLL[i][j][l];
        }
    }

    FOR (i, j)
    {
        amrex::Real ricci_hat = 0;
        FOR (k)
        {
            // We call this ricci_hat rather than ricci_tilde as we have
            // replaced what should be \tilde{Gamma} with \hat{Gamma} in
            // order to avoid adding terms that cancel later on
            ricci_hat += 0.5 * (vars.h(k, i) * d1.Gamma(k, j) +
                                vars.h(k, j) * d1.Gamma(k, i));
            ricci_hat += 0.5 * vars.Gamma(k) * d1.h(i, j, k);
            FOR (l)
            {
                ricci_hat += -0.5 * h_UU[k][l] * d2_h[i][j][k][l] +
                             (chris.ULL[k][l][i] * chris_LLU[j][k][l] +
                              chris.ULL[k][l][j] * chris_LLU[i][k][l] +
                              chris.ULL[k][i][l] * chris_LLU[k][j][l]);
            }
        }

        amrex::Real ricci_chi =
            0.5 * ((GR_SPACEDIM - 2) * covdtilde2chi[i][j] +
                   vars.h(i, j) * boxtildechi -
                   ((GR_SPACEDIM - 2) * d1.chi(i) * d1.chi(j) +
                    GR_SPACEDIM * vars.h(i, j) * dchi_dot_dchi) /
                       (2 * vars.chi()));

        amrex::Real z_terms = compute_z_terms(i, j, Z_over_chi, vars, d1.chi());

        out.LL[i][j] =
            (ricci_chi + vars.chi() * ricci_hat + z_terms) / vars.chi();
    }

    out.scalar = vars.chi() * TensorAlgebra::compute_trace(out.LL, h_UU);

    return out;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor<2, amrex::Real>
compute_d1_chris_contracted(const Tensor<2, amrex::Real> &h_UU,
                            const CCZ4D1Vars &d1,
                            const Tensor<4, amrex::Real> &d2_h)
{
    Tensor<2, amrex::Real> d1_chris_contracted = 0.0;
    FOR (i, j)
    {
        FOR (m, n, p)
        {
            amrex::Real d1_terms = 0.0;
            FOR (q, r)
            {
                d1_terms += -h_UU[q][r] * (d1.h(n, q, j) * d1.h(m, p, r) +
                                           d1.h(m, n, j) * d1.h(p, q, r));
            }
            d1_chris_contracted[i][j] +=
                h_UU[i][m] * h_UU[n][p] * (d2_h[m][n][j][p] + d1_terms);
        }
    }
    return d1_chris_contracted;
}

// This function allows adding arbitrary multiples of D_{(i}Z_{j)}
// to the Ricci scalar rather than the default of 2 in compute_ricci_Z
AMREX_GPU_DEVICE AMREX_FORCE_INLINE ricci_t compute_ricci_Z_general(
    const CCZ4Vars &vars, const CCZ4D1Vars &d1,
    const Tensor<2, amrex::Real> &d2_chi, const Tensor<4, amrex::Real> &d2_h,
    const Tensor<2, amrex::Real> &h_UU, const chris_t &chris,
    const double dZ_coeff)
{
    // get contributions from conformal metric and factor with zero Z vector
    Tensor<1, amrex::Real> zero_Z = 0.;
    auto ricci = compute_ricci_Z(vars, d1, d2_chi, d2_h, h_UU, chris, zero_Z);

    // need to add term to correct for d1.Gamma (includes Z contribution)
    // and Gamma in ricci_hat
    auto d1_chris_contracted = compute_d1_chris_contracted(h_UU, d1, d2_h);
    Tensor<1, amrex::Real> Z_over_chi;
    FOR (i)
    {
        Z_over_chi[i] = 0.5 * (vars.Gamma(i) - chris.contracted[i]);
    }
    FOR (i, j)
    {
        FOR (m)
        {
            // This corrects for the \hat{Gamma}s in ricci_hat
            ricci.LL[i][j] +=
                (1. - 0.5 * dZ_coeff) * 0.5 *
                (vars.h(m, i) * (d1_chris_contracted[m][j] - d1.Gamma(m, j)) +
                 vars.h(m, j) * (d1_chris_contracted[m][i] - d1.Gamma(m, i)) +
                 (chris.contracted[m] - vars.Gamma(m)) * d1.h(i, j, m));
        }
        amrex::Real z_terms = compute_z_terms(i, j, Z_over_chi, vars, d1.chi());
        ricci.LL[i][j] += 0.5 * dZ_coeff * z_terms / vars.chi();
    }
    ricci.scalar = vars.chi() * TensorAlgebra::compute_trace(ricci.LL, h_UU);
    return ricci;
}

// This function returns the pure Ricci scalar with no contribution from the
// Z vector - used e.g. in the constraint calculations.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE ricci_t compute_ricci(
    const CCZ4Vars &vars, const CCZ4D1Vars &d1,
    const Tensor<2, amrex::Real> &d2_chi, const Tensor<4, amrex::Real> &d2_h,
    const Tensor<2, amrex::Real> &h_UU, const chris_t &chris)
{
    return compute_ricci_Z_general(vars, d1, d2_chi, d2_h, h_UU, chris, 0.0);
}

} // namespace CCZ4Geometry

#endif /* CCZ4GEOMETRY_HPP_ */
