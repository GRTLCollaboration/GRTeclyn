/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// This file calculates CCZ4 geometric quantities (or a similar 3+1 split).
#ifndef CCZ4GEOMETRY_HPP_
#define CCZ4GEOMETRY_HPP_

#include "CCZ4Vars.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "TensorAlgebra.hpp"

//! A structure for the decomposed elements of the Energy Momentum Tensor in
//! 3+1D
struct emtensor_t
{
    Tensor::Rank2 S{}; //!< S_ij = T_ij
    Tensor::Rank1 j{}; //!< j_i = T_ia_n^a
    amrex::Real trS{}; //!< trS = S^i_i
    amrex::Real rho{}; //!< rho = T_ab n^a n^b
};

struct ricci_t
{
    Tensor::Rank2 LL{};   // Ricci with two indices down
    amrex::Real scalar{}; // Ricci scalar
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

AMREX_GPU_DEVICE AMREX_FORCE_INLINE TensorArray::Rank2
compute_inverse_metric(const CCZ4Vars &vars)
{
    amrex::Real det_h         = compute_metric_determinant(vars);
    amrex::Real det_h_inverse = 1. / det_h;
    TensorArray::Rank2 h_UU{};

    h_UU(0, 0) = (vars.h(1, 1) * vars.h(2, 2) - vars.h(1, 2) * vars.h(2, 1)) *
                 det_h_inverse;
    h_UU(0, 1) = (vars.h(2, 0) * vars.h(1, 2) - vars.h(1, 0) * vars.h(2, 2)) *
                 det_h_inverse;
    h_UU(0, 2) = (vars.h(1, 0) * vars.h(2, 1) - vars.h(2, 0) * vars.h(1, 1)) *
                 det_h_inverse;
    h_UU(1, 1) = (vars.h(0, 0) * vars.h(2, 2) - vars.h(2, 0) * vars.h(0, 2)) *
                 det_h_inverse;
    h_UU(1, 2) = (vars.h(0, 1) * vars.h(2, 0) - vars.h(0, 0) * vars.h(2, 1)) *
                 det_h_inverse;
    h_UU(2, 2) = (vars.h(0, 0) * vars.h(1, 1) - vars.h(0, 1) * vars.h(1, 0)) *
                 det_h_inverse;
    h_UU(1, 0) = h_UU(0, 1);
    h_UU(2, 0) = h_UU(0, 2);
    h_UU(2, 1) = h_UU(1, 2);

    return h_UU;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor::Rank2
compute_inverse_metric_test(const CCZ4Vars &vars)
{
    amrex::Real det_h         = compute_metric_determinant(vars);
    amrex::Real det_h_inverse = 1. / det_h;
    Tensor::Rank2 h_UU{};

    h_UU(0, 0) = (vars.h(1, 1) * vars.h(2, 2) - vars.h(1, 2) * vars.h(2, 1)) *
                 det_h_inverse;
    h_UU(0, 1) = (vars.h(2, 0) * vars.h(1, 2) - vars.h(1, 0) * vars.h(2, 2)) *
                 det_h_inverse;
    h_UU(0, 2) = (vars.h(1, 0) * vars.h(2, 1) - vars.h(2, 0) * vars.h(1, 1)) *
                 det_h_inverse;
    h_UU(1, 1) = (vars.h(0, 0) * vars.h(2, 2) - vars.h(2, 0) * vars.h(0, 2)) *
                 det_h_inverse;
    h_UU(1, 2) = (vars.h(0, 1) * vars.h(2, 0) - vars.h(0, 0) * vars.h(2, 1)) *
                 det_h_inverse;
    h_UU(2, 2) = (vars.h(0, 0) * vars.h(1, 1) - vars.h(0, 1) * vars.h(1, 0)) *
                 det_h_inverse;
    h_UU(1, 0) = h_UU(0, 1);
    h_UU(2, 0) = h_UU(0, 2);
    h_UU(2, 1) = h_UU(1, 2);

    return h_UU;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_trace_A(const CCZ4Vars &vars)
{
    TensorArray::Rank2 inverse_metric = compute_inverse_metric(vars);
    amrex::Real trace_A               = 0.0;
    FOR (i, j)
    {
        trace_A += inverse_metric(i, j) * vars.A(i, j);
    }
    return trace_A;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE TensorArray::Rank2
compute_A_UU(const CCZ4Vars &vars, const TensorArray::Rank2 &inverse_metric)
{
    TensorArray::Rank2 A_UU{};
    FOR (i, j)
    {
        A_UU(i, j) = 0.0;
        FOR (k, l)
        {
            A_UU(i, j) +=
                inverse_metric(i, k) * inverse_metric(j, l) * vars.A(k, l);
        }
    }
    return A_UU;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor::Rank2
compute_A_UU(const CCZ4Vars &vars, const Tensor::Rank2 &inverse_metric)
{
    Tensor::Rank2 A_UU{};
    FOR (i, j)
    {
        A_UU(i, j) = 0.0;
        FOR (k, l)
        {
            A_UU(i, j) +=
                inverse_metric(i, k) * inverse_metric(j, l) * vars.A(k, l);
        }
    }
    return A_UU;
}

// This is A_ij A^ij
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_Aij_squared(const CCZ4Vars &vars)
{
    TensorArray::Rank2 inverse_metric = compute_inverse_metric(vars);
    amrex::Real Aij_squared           = 0.0;
    FOR (i, j, k, l)
    {
        Aij_squared += inverse_metric(i, k) * inverse_metric(j, l) *
                       vars.A(i, j) * vars.A(k, l);
    }
    return Aij_squared;
}

// This is A_ij A^ij
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_Aij_squared(const CCZ4Vars &vars,
                    const TensorArray::Rank2 &inverse_metric)
{
    amrex::Real Aij_squared = 0.0;
    FOR (i, j, k, l)
    {
        Aij_squared += inverse_metric(i, k) * inverse_metric(j, l) *
                       vars.A(i, j) * vars.A(k, l);
    }
    return Aij_squared;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_Aij_squared(const CCZ4Vars &vars, const Tensor::Rank2 &inverse_metric)
{
    amrex::Real Aij_squared = 0.0;
    FOR (i, j, k, l)
    {
        Aij_squared += inverse_metric(i, k) * inverse_metric(j, l) *
                       vars.A(i, j) * vars.A(k, l);
    }
    return Aij_squared;
}

// This is A_ij A^ij
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_Aij_squared_with_A_UU(const CCZ4Vars &vars,
                              const TensorArray::Rank2 &A_UU)
{
    amrex::Real Aij_squared = 0.0;
    FOR (i, j)
    {
        Aij_squared += A_UU(i, j) * vars.A(i, j);
    }
    return Aij_squared;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_Aij_squared_with_A_UU(const CCZ4Vars &vars, const Tensor::Rank2 &A_UU)
{
    amrex::Real Aij_squared = 0.0;
    FOR (i, j)
    {
        Aij_squared += A_UU(i, j) * vars.A(i, j);
    }
    return Aij_squared;
}

/// Computes the conformal christoffel symbol - using AMReX Arrays

// AMREX_GPU_DEVICE AMREX_FORCE_INLINE chris_t compute_christoffel(
//     const amrex::Array2D<amrex::Real, 0, NUM_SYM_IDXS - 1, 0, AMREX_SPACEDIM
//     - 1>
//         &d1_h,
//     const TensorArray::Rank1Sym &h_UU)
// {
//     chris_t out{};

//     FOR (i, j, k)
//     {
//         out.LLL(i, j, k) =
//             0.5 * (d1_h(VAR_IDX0(j, i), k) + d1_h(VAR_IDX0(k, i), j) -
//                    d1_h(VAR_IDX0(j, k), i));
//     }

//     // FOR (i)
//     // {
//     //     out.LLL(i, i, i) = 0.5 * d1_h(SYMM_IDX(i, i), i);
//     // }

//     // out.LLL(0, 0, 1) = 0.5 * d1_h(0, 1);
//     // out.LLL(0, 0, 2) = 0.5 * d1_h(0, 2);

//     // out.LLL(0, 1, 0) = 0.5 * d1_h(0, 1);
//     // out.LLL(0, 2, 0) = 0.5 * d1_h(0, 2);

//     // out.LLL(0, 1, 1) = d1_h(1, 1) - 0.5 * d1_h(4, 0);
//     // out.LLL(0, 2, 2) = d1_h(2, 2) - 0.5 * d1_h(5, 0);

//     // out.LLL(1, 0, 0) = d1_h(1, 0) - 0.5 * d1_h(0, 1);
//     // out.LLL(1, 2, 2) = d1_h(4, 2) - 0.5 * d1_h(5, 1);

//     // out.LLL(1, 0, 1) = 0.5 * d1_h(3, 0);
//     // out.LLL(1, 1, 0) = 0.5 * d1_h(3, 0);

//     // out.LLL(1, 1, 2) = 0.5 * d1_h(3, 2);
//     // out.LLL(1, 2, 1) = 0.5 * d1_h(3, 2);

//     // out.LLL(2, 0, 0) = d1_h(2, 0) - 0.5 * d1_h(0, 2);
//     // out.LLL(2, 1, 1) = d1_h(4, 1) - 0.5 * d1_h(3, 2);

//     // out.LLL(2, 1, 2) = 0.5 * d1_h(5, 1);
//     // out.LLL(2, 2, 1) = 0.5 * d1_h(5, 1);

//     // out.LLL(2, 2, 0) = 0.5 * d1_h(5, 0);
//     // out.LLL(2, 0, 2) = 0.5 * d1_h(5, 0);

//     // ////These have all different indices

//     // out.LLL(0, 1, 2) = 0.5 * (d1_h(1, 2) + d1_h(2, 0) + d1_h(4, 0));
//     // out.LLL(0, 2, 1) = out.LLL(0, 1, 2);

//     // out.LLL(1, 0, 2) = 0.5 * (d1_h(1, 2) + d1_h(4, 0) - d1_h(2, 0));
//     // out.LLL(1, 2, 0) = out.LLL(1, 0, 2);

//     // out.LLL(2, 0, 1) = 0.5 * (d1_h(3, 1) + d1_h(4, 0) - d1_h(1, 2));
//     // out.LLL(2, 1, 0) = out.LLL(2, 0, 1);

//     FOR (i, j)
//     {
//         out.ULL(0, i, j) = h_UU(0) * out.LLL(0, i, j) +
//                            h_UU(1) * out.LLL(1, i, j) +
//                            h_UU(2) * out.LLL(2, i, j);
//         out.ULL(1, i, j) = h_UU(1) * out.LLL(0, i, j) +
//                            h_UU(3) * out.LLL(1, i, j) +
//                            h_UU(4) * out.LLL(2, i, j);
//         out.ULL(2, i, j) = h_UU(2) * out.LLL(0, i, j) +
//                            h_UU(4) * out.LLL(1, i, j) +
//                            h_UU(5) * out.LLL(2, i, j);
//     }

//     // FOR (i, j, k)
//     // {
//     //     out.ULL(i, j, k) = 0;
//     //     FOR (l)
//     //     {
//     //         out.ULL(i, j, k) += h_UU(SYMM_IDX(i, l)) * out.LLL(l, j, k);
//     //     }
//     // }
//     FOR (i)
//     {
//         out.contracted(i) =
//             h_UU(0) * out.ULL(i, 0, 0) + h_UU(1) * out.ULL(i, 0, 1) +
//             h_UU(2) * out.ULL(i, 0, 2) + h_UU(1) * out.ULL(i, 1, 0) +
//             h_UU(3) * out.ULL(i, 1, 1) + h_UU(4) * out.ULL(i, 1, 2) +
//             h_UU(2) * out.ULL(i, 2, 0) + h_UU(4) * out.ULL(i, 2, 1) +
//             h_UU(5) * out.ULL(i, 2, 2);
//     }

//     return out;
// }

/// Computes the conformal christoffel symbol - using tensors
AMREX_GPU_DEVICE AMREX_FORCE_INLINE chris_t
compute_christoffel(const amrex::Array2D<amrex::Real, 0, NUM_SYM_IDXS - 1, 0,
                                         AMREX_SPACEDIM - 1> &d1_h,
                    const TensorArray::Rank2 &h_UU)
{
    chris_t out{};

    FOR (i, j, k)
    {
        out.LLL(i, j, k) =
            0.5 * (d1_h(VAR_IDX0(j, i), k) + d1_h(VAR_IDX0(k, i), j) -
                   d1_h(VAR_IDX0(j, k), i));
    }

    FOR (i, j, k)
    {
        out.ULL(i, j, k) = 0;
        FOR (l)
        {
            out.ULL(i, j, k) += h_UU(i, l) * out.LLL(l, j, k);
        }
    }
    FOR (i)
    {
        out.contracted(i) = 0;
        FOR (j, k)
        {
            out.contracted(i) += h_UU(j, k) * out.ULL(i, j, k);
        }
    }

    return out;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE chris_t
compute_christoffel(const Tensor::Sym12Rank3 &d1_h, const Tensor::Rank2 &h_UU)
{
    chris_t out{};

    FOR (i, j, k)
    {
        out.LLL(i, j, k) =
            0.5 * (d1_h(j, i, k) + d1_h(k, i, j) - d1_h(j, k, i));
    }

    FOR (i, j, k)
    {
        out.ULL(i, j, k) = 0;
        FOR (l)
        {
            out.ULL(i, j, k) += h_UU(i, l) * out.LLL(l, j, k);
        }
    }
    FOR (i)
    {
        out.contracted(i) = 0;
        FOR (j, k)
        {
            out.contracted(i) += h_UU(j, k) * out.ULL(i, j, k);
        }
    }

    return out;
}

/// Computes the conformal christoffel symbol
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE TensorArray::Rank3
compute_phys_chris(const CCZ4Vars &vars, const TensorArray::Rank1 &d1_chi,
                   const TensorArray::Rank2 &h_UU,
                   const TensorArray::Rank3 &chris_ULL)
{
    using namespace TensorAlgebra;
    TensorArray::Rank3 chris_phys{};
    FOR (i, j, k)
    {
        chris_phys(i, j, k) =
            chris_ULL(i, j, k) -
            0.5 / vars.chi() *
                (delta(i, k) * d1_chi(j) + delta(i, j) * d1_chi(k));
        FOR (m)
        {
            chris_phys(i, j, k) +=
                0.5 / vars.chi() * vars.h(j, k) * h_UU(i, m) * d1_chi(m);
        }
    }
    return chris_phys;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE Tensor::Rank3
compute_phys_chris(const CCZ4Vars &vars, const Tensor::Rank1 &d1_chi,
                   const Tensor::Rank2 &h_UU, const Tensor::Rank3 &chris_ULL)
{
    using namespace TensorAlgebra;
    Tensor::Rank3 chris_phys{};
    FOR (i, j, k)
    {
        chris_phys(i, j, k) =
            chris_ULL(i, j, k) -
            0.5 / vars.chi() *
                (delta(i, k) * d1_chi(j) + delta(i, j) * d1_chi(k));
        FOR (m)
        {
            chris_phys(i, j, k) +=
                0.5 / vars.chi() * vars.h(j, k) * h_UU(i, m) * d1_chi(m);
        }
    }
    return chris_phys;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_divshift(TensorArray::Rank2 &d1_shift)
{
    amrex::Real divshift = 0.;
    FOR (i)
        divshift += d1_shift(i, i);
    return divshift;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_divshift(Tensor::Rank2 &d1_shift)
{
    amrex::Real divshift = 0.;
    FOR (i)
        divshift += d1_shift(i, i);
    return divshift;
}

/// Removes the trace of a 2-Tensor with lower indices given a metric and an
/// inverse metric.

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
make_trace_free(TensorArray::Rank2 &tensor_LL, const CCZ4Vars vars,
                const TensorArray::Rank2 &inverse_metric)
{
    auto trace = TensorAlgebra::compute_trace(tensor_LL, inverse_metric);
    double one_over_gr_spacedim = 1. / ((double)GR_SPACEDIM);
    FOR (i, j)
    {
        tensor_LL(i, j) -= one_over_gr_spacedim * vars.h(i, j) * trace;
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
make_trace_free(TensorArray::Rank2 &tensor_LL, const CCZ4Vars vars,
                const Tensor::Rank2 &inverse_metric)
{
    auto trace = TensorAlgebra::compute_trace(tensor_LL, inverse_metric);
    double one_over_gr_spacedim = 1. / ((double)GR_SPACEDIM);
    FOR (i, j)
    {
        tensor_LL(i, j) -= one_over_gr_spacedim * vars.h(i, j) * trace;
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
make_trace_free(Tensor::Rank2 &tensor_LL, const CCZ4Vars vars,
                const Tensor::Rank2 &inverse_metric)
{
    auto trace = TensorAlgebra::compute_trace(tensor_LL, inverse_metric);
    double one_over_gr_spacedim = 1. / ((double)GR_SPACEDIM);
    FOR (i, j)
    {
        tensor_LL(i, j) -= one_over_gr_spacedim * vars.h(i, j) * trace;
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
make_trace_free(Tensor::Rank2 &tensor_LL, const CCZ4Vars vars,
                const TensorArray::Rank2 &inverse_metric)
{
    auto trace = TensorAlgebra::compute_trace(tensor_LL, inverse_metric);
    double one_over_gr_spacedim = 1. / ((double)GR_SPACEDIM);
    FOR (i, j)
    {
        tensor_LL(i, j) -= one_over_gr_spacedim * vars.h(i, j) * trace;
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_z_terms(const int i, const int j, const TensorArray::Rank1 &Z_over_chi,
                const CCZ4Vars &vars, const TensorArray::Rank1 &d1_chi)
{
    amrex::Real out = 0.;
    FOR (k)
    {
        out += Z_over_chi(k) *
               (vars.h(i, k) * d1_chi(j) + vars.h(j, k) * d1_chi(i) -
                vars.h(i, j) * d1_chi(k));
    }
    return out;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_z_terms(const int i, const int j, const Tensor::Rank1 &Z_over_chi,
                const CCZ4Vars &vars, const Tensor::Rank1 &d1_chi)
{
    amrex::Real out = 0.;
    FOR (k)
    {
        out += Z_over_chi(k) *
               (vars.h(i, k) * d1_chi(j) + vars.h(j, k) * d1_chi(i) -
                vars.h(i, j) * d1_chi(k));
    }
    return out;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_ricci_hat(const int i, const int j, const CCZ4Vars &vars,
                  const TensorArray::Rank2 &d1_Gamma,
                  const amrex::Array2D<amrex::Real, 0, NUM_SYM_IDXS - 1, 0,
                                       AMREX_SPACEDIM - 1> &d1_h,
                  const TensorArray::Rank2Sym &d2_h,
                  const TensorArray::Rank2 &h_UU, const chris_t &chris)
{

    amrex::Real ricci_hat = 0;
    int idx1              = VAR_IDX0(i, j);

    FOR (k)
    {

        // We call this ricci_hat rather than ricci_tilde as we have
        // replaced what should be \tilde{Gamma} with \hat{Gamma} in
        // order to avoid adding terms that cancel later on
        ricci_hat += 0.5 * (vars.h(k, i) * d1_Gamma(k, j) +
                            vars.h(k, j) * d1_Gamma(k, i));
        ricci_hat += 0.5 * vars.Gamma(k) * d1_h(idx1, k);

        FOR (l)
        {
            amrex::Real chris_LLU_jkl = 0.0;
            amrex::Real chris_LLU_ikl = 0.0;
            amrex::Real chris_LLU_kjl = 0.0;

            FOR (m)
            {
                // jkl
                chris_LLU_jkl += h_UU(l, m) * chris.LLL(j, k, m);

                // ikl
                chris_LLU_ikl += h_UU(l, m) * chris.LLL(i, k, m);

                // kjl
                chris_LLU_kjl += h_UU(l, m) * chris.LLL(k, j, m);
            }

            int idx2 = VAR_IDX0(k, l);

            ricci_hat += -0.5 * h_UU(k, l) * d2_h(idx1, idx2) +
                         chris.ULL(k, l, i) * chris_LLU_jkl +
                         chris.ULL(k, l, j) * chris_LLU_ikl +
                         chris.ULL(k, i, l) * chris_LLU_kjl;
        }
    }

    return ricci_hat;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_ricci_hat(const int i, const int j, const CCZ4Vars &vars,
                  const Tensor::Rank2 &d1_Gamma, const Tensor::Sym12Rank3 &d1_h,
                  const Tensor::Sym12Sym34Rank4 &d2_h,
                  const Tensor::Rank2 &h_UU, const chris_t &chris)
{

    amrex::Real ricci_hat = 0;
    int idx1              = VAR_IDX0(i, j);

    FOR (k)
    {

        // We call this ricci_hat rather than ricci_tilde as we have
        // replaced what should be \tilde{Gamma} with \hat{Gamma} in
        // order to avoid adding terms that cancel later on
        ricci_hat += 0.5 * (vars.h(k, i) * d1_Gamma(k, j) +
                            vars.h(k, j) * d1_Gamma(k, i));
        ricci_hat += 0.5 * vars.Gamma(k) * d1_h(idx1, k);

        FOR (l)
        {
            amrex::Real chris_LLU_jkl = 0.0;
            amrex::Real chris_LLU_ikl = 0.0;
            amrex::Real chris_LLU_kjl = 0.0;

            FOR (m)
            {
                // jkl
                chris_LLU_jkl += h_UU(l, m) * chris.LLL(j, k, m);

                // ikl
                chris_LLU_ikl += h_UU(l, m) * chris.LLL(i, k, m);

                // kjl
                chris_LLU_kjl += h_UU(l, m) * chris.LLL(k, j, m);
            }

            int idx2 = VAR_IDX0(k, l);

            ricci_hat += -0.5 * h_UU(k, l) * d2_h(idx1, idx2) +
                         chris.ULL(k, l, i) * chris_LLU_jkl +
                         chris.ULL(k, l, j) * chris_LLU_ikl +
                         chris.ULL(k, i, l) * chris_LLU_kjl;
        }
    }

    return ricci_hat;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE ricci_t compute_ricci_Z(
    const CCZ4Vars &vars, const TensorArray::Rank1 &d1_chi,
    const TensorArray::Rank2 &d1_Gamma,
    const amrex::Array2D<amrex::Real, 0, NUM_SYM_IDXS - 1, 0,
                         AMREX_SPACEDIM - 1> &d1_h,
    const TensorArray::Rank2Sym &d2_h, const TensorArray::Rank1Sym &d2_chi,
    const TensorArray::Rank2 &h_UU, const chris_t &chris,
    const TensorArray::Rank1 &Z_over_chi)
{
    ricci_t out;

    TensorArray::Rank2 covdtilde2chi{};
    FOR (k, l)
    {
        covdtilde2chi(k, l) = d2_chi(VAR_IDX0(k, l));
        FOR (m)
        {
            covdtilde2chi(k, l) -= chris.ULL(m, k, l) * d1_chi(m);
        }
    }

    amrex::Real boxtildechi   = 0.;
    amrex::Real dchi_dot_dchi = 0.;
    FOR (i, j)
    {
        boxtildechi   += covdtilde2chi(i, j) * h_UU(i, j);
        dchi_dot_dchi += d1_chi(i) * d1_chi(j) * h_UU(i, j);
    }

    FOR (i, j)
    {

        amrex::Real ricci_hat =
            compute_ricci_hat(i, j, vars, d1_Gamma, d1_h, d2_h, h_UU, chris);

        amrex::Real ricci_chi =
            0.5 * ((GR_SPACEDIM - 2) * covdtilde2chi(i, j) +
                   vars.h(i, j) * boxtildechi -
                   ((GR_SPACEDIM - 2) * d1_chi(i) * d1_chi(j) +
                    GR_SPACEDIM * vars.h(i, j) * dchi_dot_dchi) /
                       (2 * vars.chi()));

        amrex::Real z_terms = compute_z_terms(i, j, Z_over_chi, vars, d1_chi);

        out.LL(i, j) =
            (ricci_chi + vars.chi() * ricci_hat + z_terms) / vars.chi();
    }

    out.scalar = vars.chi() * TensorAlgebra::compute_trace(out.LL, h_UU);

    return out;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE ricci_t
compute_ricci_Z(const CCZ4Vars &vars, const Tensor::Rank1 &d1_chi,
                const Tensor::Rank2 &d1_Gamma, const Tensor::Sym12Rank3 &d1_h,
                const Tensor::Sym12Sym34Rank4 &d2_h,
                const Tensor::Sym12Rank2 &d2_chi, const Tensor::Rank2 &h_UU,
                const chris_t &chris, const Tensor::Rank1 &Z_over_chi)
{
    ricci_t out;

    Tensor::Rank2 covdtilde2chi{};
    FOR (k, l)
    {
        covdtilde2chi(k, l) = d2_chi(k, l);
        FOR (m)
        {
            covdtilde2chi(k, l) -= chris.ULL(m, k, l) * d1_chi(m);
        }
    }

    amrex::Real boxtildechi   = 0.;
    amrex::Real dchi_dot_dchi = 0.;
    FOR (i, j)
    {
        boxtildechi   += covdtilde2chi(i, j) * h_UU(i, j);
        dchi_dot_dchi += d1_chi(i) * d1_chi(j) * h_UU(i, j);
    }

    FOR (i, j)
    {

        amrex::Real ricci_hat =
            compute_ricci_hat(i, j, vars, d1_Gamma, d1_h, d2_h, h_UU, chris);

        amrex::Real ricci_chi =
            0.5 * ((GR_SPACEDIM - 2) * covdtilde2chi(i, j) +
                   vars.h(i, j) * boxtildechi -
                   ((GR_SPACEDIM - 2) * d1_chi(i) * d1_chi(j) +
                    GR_SPACEDIM * vars.h(i, j) * dchi_dot_dchi) /
                       (2 * vars.chi()));

        amrex::Real z_terms = compute_z_terms(i, j, Z_over_chi, vars, d1_chi);

        out.LL(i, j) =
            (ricci_chi + vars.chi() * ricci_hat + z_terms) / vars.chi();
    }

    out.scalar = vars.chi() * TensorAlgebra::compute_trace(out.LL, h_UU);

    return out;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE TensorArray::Rank2
compute_d1_chris_contracted(
    const TensorArray::Rank2 &h_UU,
    const amrex::Array2D<amrex::Real, 0, NUM_SYM_IDXS - 1, 0,
                         AMREX_SPACEDIM - 1> &d1_h,
    const TensorArray::Rank2Sym &d2_h)
{
    TensorArray::Rank2 d1_chris_contracted{};

    FOR (i, j)
    {
        d1_chris_contracted(i, j) = 0.;

        FOR (m, n, p)
        {
            amrex::Real d1_terms = 0.0;
            FOR (q, r)
            {
                d1_terms += -h_UU(q, r) *
                            (d1_h(VAR_IDX0(n, q), j) * d1_h(VAR_IDX0(m, p), r) +
                             d1_h(VAR_IDX0(m, n), j) * d1_h(VAR_IDX0(p, q), r));
            }

            d1_chris_contracted(i, j) +=
                h_UU(i, m) * h_UU(n, p) *
                (d2_h(VAR_IDX0(m, n), VAR_IDX0(j, p)) + d1_terms);
        }
    }
    return d1_chris_contracted;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor::Rank2
compute_d1_chris_contracted(const Tensor::Rank2 &h_UU,
                            const Tensor::Sym12Rank3 &d1_h,
                            const Tensor::Sym12Sym34Rank4 &d2_h)
{
    Tensor::Rank2 d1_chris_contracted{};

    FOR (i, j)
    {
        d1_chris_contracted(i, j) = 0.;

        FOR (m, n, p)
        {
            amrex::Real d1_terms = 0.0;
            FOR (q, r)
            {
                d1_terms += -h_UU(q, r) * (d1_h(n, q, j) * d1_h(m, p, r) +
                                           d1_h(m, n, j) * d1_h(p, q, r));
            }

            d1_chris_contracted(i, j) +=
                h_UU(i, m) * h_UU(n, p) * (d2_h(m, n, j, p) + d1_terms);
        }
    }
    return d1_chris_contracted;
}

// This function allows adding arbitrary multiples of D_{(i}Z_{j)}
// to the Ricci scalar rather than the default of 2 in compute_ricci_Z
AMREX_GPU_DEVICE AMREX_FORCE_INLINE ricci_t compute_ricci_Z_general(
    const CCZ4Vars &vars, const TensorArray::Rank1 &d1_chi,
    const TensorArray::Rank2 &d1_Gamma,
    const amrex::Array2D<amrex::Real, 0, NUM_SYM_IDXS - 1, 0,
                         AMREX_SPACEDIM - 1> &d1_h,
    const TensorArray::Rank1Sym &d2_chi, const TensorArray::Rank2Sym &d2_h,
    const TensorArray::Rank2 &h_UU, const chris_t &chris, const double dZ_coeff)
{
    // get contributions from conformal metric and factor with zero Z vector
    TensorArray::Rank1 zero_Z{};

    FOR (i)
    {
        zero_Z(i) = 0.;
    }

    auto ricci = compute_ricci_Z(vars, d1_chi, d1_Gamma, d1_h, d2_h, d2_chi,
                                 h_UU, chris, zero_Z);

    // need to add term to correct for d1.Gamma (includes Z contribution)
    // and Gamma in ricci_hat
    auto d1_chris_contracted = compute_d1_chris_contracted(h_UU, d1_h, d2_h);
    TensorArray::Rank1 Z_over_chi{};
    FOR (i)
    {
        Z_over_chi(i) = 0.5 * (vars.Gamma(i) - chris.contracted(i));
    }
    FOR (i, j)
    {
        FOR (m)
        {
            // This corrects for the \hat{Gamma}s in ricci_hat
            ricci.LL(i, j) +=
                (1. - 0.5 * dZ_coeff) * 0.5 *
                (vars.h(m, i) * (d1_chris_contracted(m, j) - d1_Gamma(m, j)) +
                 vars.h(m, j) * (d1_chris_contracted(m, i) - d1_Gamma(m, i)) +
                 (chris.contracted(m) - vars.Gamma(m)) *
                     d1_h(VAR_IDX0(i, j), m));
        }
        amrex::Real z_terms  = compute_z_terms(i, j, Z_over_chi, vars, d1_chi);
        ricci.LL(i, j)      += 0.5 * dZ_coeff * z_terms / vars.chi();
    }
    ricci.scalar = vars.chi() * TensorAlgebra::compute_trace(ricci.LL, h_UU);
    return ricci;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE ricci_t compute_ricci_Z_general(
    const CCZ4Vars &vars, const Tensor::Rank1 &d1_chi,
    const Tensor::Rank2 &d1_Gamma, const Tensor::Sym12Rank3 &d1_h,
    const Tensor::Sym12Rank2 &d2_chi, const Tensor::Sym12Sym34Rank4 &d2_h,
    const Tensor::Rank2 &h_UU, const chris_t &chris, const double dZ_coeff)
{
    // get contributions from conformal metric and factor with zero Z vector
    Tensor::Rank1 zero_Z{};

    FOR (i)
    {
        zero_Z(i) = 0.;
    }

    auto ricci = compute_ricci_Z(vars, d1_chi, d1_Gamma, d1_h, d2_h, d2_chi,
                                 h_UU, chris, zero_Z);

    // need to add term to correct for d1.Gamma (includes Z contribution)
    // and Gamma in ricci_hat
    auto d1_chris_contracted = compute_d1_chris_contracted(h_UU, d1_h, d2_h);
    Tensor::Rank1 Z_over_chi{};
    FOR (i)
    {
        Z_over_chi(i) = 0.5 * (vars.Gamma(i) - chris.contracted(i));
    }
    FOR (i, j)
    {
        FOR (m)
        {
            // This corrects for the \hat{Gamma}s in ricci_hat
            ricci.LL(i, j) +=
                (1. - 0.5 * dZ_coeff) * 0.5 *
                (vars.h(m, i) * (d1_chris_contracted(m, j) - d1_Gamma(m, j)) +
                 vars.h(m, j) * (d1_chris_contracted(m, i) - d1_Gamma(m, i)) +
                 (chris.contracted(m) - vars.Gamma(m)) * d1_h(i, j, m));
        }
        amrex::Real z_terms  = compute_z_terms(i, j, Z_over_chi, vars, d1_chi);
        ricci.LL(i, j)      += 0.5 * dZ_coeff * z_terms / vars.chi();
    }
    ricci.scalar = vars.chi() * TensorAlgebra::compute_trace(ricci.LL, h_UU);
    return ricci;
}

// This function returns the pure Ricci scalar with no contribution from the
// Z vector - used e.g. in the constraint calculations.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE ricci_t compute_ricci(
    const CCZ4Vars &vars, const TensorArray::Rank1 &d1_chi,
    const TensorArray::Rank2 &d1_Gamma,
    const amrex::Array2D<amrex::Real, 0, NUM_SYM_IDXS - 1, 0,
                         AMREX_SPACEDIM - 1> &d1_h,
    const TensorArray::Rank1Sym &d2_chi, const TensorArray::Rank2Sym &d2_h,
    const TensorArray::Rank2 &h_UU, const chris_t &chris)
{
    return compute_ricci_Z_general(vars, d1_chi, d1_Gamma, d1_h, d2_chi, d2_h,
                                   h_UU, chris, 0.0);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE ricci_t compute_ricci(
    const CCZ4Vars &vars, const Tensor::Rank1 &d1_chi,
    const Tensor::Rank2 &d1_Gamma, const Tensor::Sym12Rank3 &d1_h,
    const Tensor::Sym12Rank2 &d2_chi, const Tensor::Sym12Sym34Rank4 &d2_h,
    const Tensor::Rank2 &h_UU, const chris_t &chris)
{
    return compute_ricci_Z_general(vars, d1_chi, d1_Gamma, d1_h, d2_chi, d2_h,
                                   h_UU, chris, 0.0);
}

} // namespace CCZ4Geometry

#endif /* CCZ4GEOMETRY_HPP_ */
