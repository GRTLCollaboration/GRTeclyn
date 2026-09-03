/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef TENSORALGEBRA_HPP_
#define TENSORALGEBRA_HPP_

// GRTeclyn includes
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"

// AMReX includes
#include <AMReX_Array.H>
#include <AMReX_REAL.H>

// System includes
#include <array>

using namespace amrex::literals;

struct chris_t
{
    Tensor::Rank3 ULL{};
    Tensor::Rank3 LLL{};
    Tensor::Rank1 contracted{};
};

namespace TensorAlgebra
{
/// Computes determinant of a symmetric 3x3 matrix
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_determinant_sym(const Tensor::Sym12Rank2 &matrix)
{
    amrex::Real det = matrix(0, 0) * matrix(1, 1) * matrix(2, 2) +
                      2 * matrix(0, 1) * matrix(0, 2) * matrix(1, 2) -
                      matrix(0, 0) * matrix(1, 2) * matrix(1, 2) -
                      matrix(1, 1) * matrix(0, 2) * matrix(0, 2) -
                      matrix(2, 2) * matrix(0, 1) * matrix(0, 1);

    return det;
}

/// Computes the determinant of a general 3x3 matrix.
/// Note: for a symmetric matrix use the simplified function
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_determinant(const Tensor::Rank2 &matrix)
{
    amrex::Real det =
        matrix(0, 0) *
            (matrix(1, 1) * matrix(2, 2) - matrix(1, 2) * matrix(2, 1)) -
        matrix(0, 1) *
            (matrix(2, 2) * matrix(1, 0) - matrix(1, 2) * matrix(2, 0)) +
        matrix(0, 2) *
            (matrix(1, 0) * matrix(2, 1) - matrix(1, 1) * matrix(2, 0));
    return det;
}

/// Computes the inverse of a symmetric 3x3 matrix directly.
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor::Sym12Rank2
compute_inverse_sym(const Tensor::Sym12Rank2 &matrix)
{
    amrex::Real deth         = compute_determinant_sym(matrix);
    amrex::Real deth_inverse = 1._rt / deth;
    Tensor::Sym12Rank2 h_UU{};
    h_UU(0, 0) = (matrix(1, 1) * matrix(2, 2) - matrix(1, 2) * matrix(1, 2)) *
                 deth_inverse;
    h_UU(0, 1) = (matrix(0, 2) * matrix(1, 2) - matrix(0, 1) * matrix(2, 2)) *
                 deth_inverse;
    h_UU(0, 2) = (matrix(0, 1) * matrix(1, 2) - matrix(0, 2) * matrix(1, 1)) *
                 deth_inverse;
    h_UU(1, 1) = (matrix(0, 0) * matrix(2, 2) - matrix(0, 2) * matrix(0, 2)) *
                 deth_inverse;
    h_UU(1, 2) = (matrix(0, 1) * matrix(0, 2) - matrix(0, 0) * matrix(1, 2)) *
                 deth_inverse;
    h_UU(2, 2) = (matrix(0, 0) * matrix(1, 1) - matrix(0, 1) * matrix(0, 1)) *
                 deth_inverse;

    return h_UU;
}

/// Computes the inverse of a general 3x3 matrix.
/// Note: for a symmetric matrix use the simplified function
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor::Rank2
compute_inverse(const Tensor::Rank2 &matrix)
{
    amrex::Real deth         = compute_determinant(matrix);
    amrex::Real deth_inverse = 1._rt / deth;
    Tensor::Rank2 h_UU{};
    h_UU(0, 0) = (matrix(1, 1) * matrix(2, 2) - matrix(1, 2) * matrix(2, 1)) *
                 deth_inverse;
    h_UU(1, 1) = (matrix(0, 0) * matrix(2, 2) - matrix(0, 2) * matrix(2, 0)) *
                 deth_inverse;
    h_UU(2, 2) = (matrix(0, 0) * matrix(1, 1) - matrix(1, 0) * matrix(0, 1)) *
                 deth_inverse;
    h_UU(1, 0) = (matrix(2, 0) * matrix(1, 2) - matrix(1, 0) * matrix(2, 2)) *
                 deth_inverse;
    h_UU(0, 1) = (matrix(0, 2) * matrix(2, 1) - matrix(0, 1) * matrix(2, 2)) *
                 deth_inverse;
    h_UU(2, 0) = (matrix(1, 0) * matrix(2, 1) - matrix(1, 1) * matrix(2, 0)) *
                 deth_inverse;
    h_UU(0, 2) = (matrix(0, 1) * matrix(1, 2) - matrix(1, 1) * matrix(0, 2)) *
                 deth_inverse;
    h_UU(2, 1) = (matrix(0, 1) * matrix(2, 0) - matrix(0, 0) * matrix(2, 1)) *
                 deth_inverse;
    h_UU(1, 2) = (matrix(1, 0) * matrix(0, 2) - matrix(0, 0) * matrix(1, 2)) *
                 deth_inverse;

    return h_UU;
}

/// Computes the trace of a 2-Tensor with lower indices given an inverse metric.
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_trace(const Tensor::Rank2 &tensor_LL,
              const Tensor::Rank2 &inverse_metric)
{
    amrex::Real trace = 0._rt;
    FOR (i, j)
    {
        trace += inverse_metric(i, j) * tensor_LL(i, j);
    }
    return trace;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_trace(const Tensor::Rank2 &tensor_LL,
              const Tensor::Sym12Rank2 &sym_inverse_metric)
{
    amrex::Real trace = 0._rt;
    FOR (i, j)
    {
        trace += sym_inverse_metric(i, j) * tensor_LL(i, j);
    }
    return trace;
}

/// Computes the trace of a 1,1 Tensor (a matrix) - no metric required.
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_trace(const Tensor::Rank2 &tensor_UL)
{
    amrex::Real trace = 0._rt;
    FOR (i)
        trace += tensor_UL(i, i);
    return trace;
}

/// Computes dot product of a vector and a covector (no metric required)
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_dot_product(const Tensor::Rank1 &vector_U,
                    const Tensor::Rank1 &covector_L)
{
    amrex::Real dot_product = 0._rt;
    FOR (i)
        dot_product += vector_U(i) * covector_L(i);
    return dot_product;
}

/// Computes dot product of two covectors given an inverse metric or
/// the dot product of two vectors given a metric.
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_dot_product(const Tensor::Rank1 &covector1_L,
                    const Tensor::Rank1 &covector2_L,
                    const Tensor::Rank2 &inverse_metric)
{
    amrex::Real dot_product = 0._rt;
    FOR (m, n)
    {
        dot_product += inverse_metric(m, n) * covector1_L(m) * covector2_L(n);
    }
    return dot_product;
}

/// Computes dot product of two covectors given an inverse metric or
/// the dot product of two vectors given a (symmetric) metric.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_dot_product(const Tensor::Rank1 &covector1_L,
                    const Tensor::Rank1 &covector2_L,
                    const Tensor::Sym12Rank2 &inverse_metric_sym)
{
    amrex::Real dot_product = 0._rt;
    FOR (m, n)
    {
        dot_product +=
            inverse_metric_sym(m, n) * covector1_L(m) * covector2_L(n);
    }
    return dot_product;
}

/// Removes the trace of a 2-Tensor with lower indices given a metric and an
/// inverse metric.  Or a Tensor with upper indices given an inverse metric and
/// a metric.

// NOLINTBEGIN
template <int size = AMREX_SPACEDIM>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE

    void
    make_trace_free(Tensor::GeneralRank<2, size, size> &tensor_LL,
                    const Tensor::GeneralRank<2, size, size> &metric,
                    const Tensor::GeneralRank<2, size, size> &inverse_metric)
// NOLINTEND

{
    auto trace                       = compute_trace(tensor_LL, inverse_metric);
    amrex::Real one_over_gr_spacedim = 1._rt / ((amrex::Real)GR_SPACEDIM);
    FOR (i, j)
    {
        tensor_LL(i, j) -= one_over_gr_spacedim * metric(i, j) * trace;
    }
}

/// Makes a 2-Tensor symmetric

// NOLINTBEGIN(readability-named-parameter)
template <int size>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
make_symmetric(Tensor::GeneralRank<2, size, size> &tensor_LL)
// NOLINTEND(readability-named-parameter)
{
    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < i; ++j)
        {
            tensor_LL(i, j) = 0.5_rt * (tensor_LL(i, j) + tensor_LL(j, i));
            tensor_LL(j, i) = tensor_LL(i, j);
        }
    }
}

/// Raises the index of a covector
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor::Rank1
raise_all(const Tensor::Rank1 &tensor_L, const Tensor::Rank2 &inverse_metric)
{
    Tensor::Rank1 tensor_U{};

    FOR (i)
    {
        tensor_U(i) = 0._rt;
    }

    FOR (i, j)
    {
        tensor_U(i) += inverse_metric(i, j) * tensor_L(j);
    }
    return tensor_U;
}

/// Raises the indices of a 2-Tensor
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor::Rank2
raise_all(const Tensor::Rank2 &tensor_LL, const Tensor::Rank2 &inverse_metric)
{
    Tensor::Rank2 tensor_UU{};

    FOR (i, j)
    {
        tensor_UU(i, j) = 0._rt;
    }

    FOR (i, j, k, l)
    {
        tensor_UU(i, j) +=
            inverse_metric(i, k) * inverse_metric(j, l) * tensor_LL(k, l);
    }
    return tensor_UU;
}

/// Lowers the indices of a vector
/// Note: same functionality as raise; included to improve readability

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor::Rank1
lower_all(const Tensor::Rank1 &tensor_U, const Tensor::Rank2 &metric)
// NOLINTEND(bugprone-easily-swappable-parameters)
{ // The code for lowering is exactly the same as for raising
    return raise_all(tensor_U, metric);
}

/// Lowers the indices of a 2-Tensor
/// Note: same functionality as raise; included to improve readability

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor::Rank2
lower_all(const Tensor::Rank2 &tensor_UU, const Tensor::Rank2 &metric)

// NOLINTEND(bugprone-easily-swappable-parameters)
{ // The code for lowering is exactly the same as for raising
    return raise_all(tensor_UU, metric);
}

/// Computes the (i,j) component of the Kronecker delta

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
constexpr int delta(int i, int j) { return static_cast<int>(i == j); }
// NOLINTEND(bugprone-easily-swappable-parameters)

[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor::Rank3 epsilon()
{
    Tensor::Rank3 epsilon(0._rt);

    epsilon(0, 1, 2) = 1.0_rt;
    epsilon(1, 2, 0) = 1.0_rt;
    epsilon(2, 0, 1) = 1.0_rt;
    epsilon(0, 2, 1) = -1.0_rt;
    epsilon(2, 1, 0) = -1.0_rt;
    epsilon(1, 0, 2) = -1.0_rt;

    return epsilon;
}

/// Computes the levi-civita symbol (4D, NB, symbol, not the Tensor)
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor::SpacetimeRank4
epsilon4D()
{
    Tensor::SpacetimeRank4 epsilon4D(0._rt);

    // Fortran order!
    epsilon4D(0, 1, 2, 3) = 1.0_rt;
    epsilon4D(0, 1, 3, 2) = -1.0_rt;
    epsilon4D(0, 3, 1, 2) = 1.0_rt;
    epsilon4D(0, 3, 2, 1) = -1.0_rt;
    epsilon4D(0, 2, 1, 3) = -1.0_rt;
    epsilon4D(0, 2, 3, 1) = 1.0_rt;

    epsilon4D(1, 0, 2, 3) = -1.0_rt;
    epsilon4D(1, 2, 0, 3) = 1.0_rt;
    epsilon4D(1, 2, 3, 0) = -1.0_rt;
    epsilon4D(1, 3, 2, 0) = 1.0_rt;
    epsilon4D(1, 3, 0, 2) = -1.0_rt;
    epsilon4D(1, 0, 3, 2) = 1.0_rt;

    epsilon4D(2, 0, 1, 3) = 1.0_rt;
    epsilon4D(2, 0, 3, 1) = -1.0_rt;
    epsilon4D(2, 3, 0, 1) = 1.0_rt;
    epsilon4D(2, 3, 1, 0) = -1.0_rt;
    epsilon4D(2, 1, 3, 0) = 1.0_rt;
    epsilon4D(2, 1, 0, 3) = -1.0_rt;

    epsilon4D(3, 0, 1, 2) = -1.0_rt;
    epsilon4D(3, 1, 0, 2) = 1.0_rt;
    epsilon4D(3, 1, 2, 0) = -1.0_rt;
    epsilon4D(3, 2, 1, 0) = 1.0_rt;
    epsilon4D(3, 2, 0, 1) = -1.0_rt;
    epsilon4D(3, 0, 2, 1) = 1.0_rt;

    return epsilon4D;
}

/// Computes the conformal christoffel symbol

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE chris_t
compute_christoffel(const Tensor::Rank3 &d1_metric, const Tensor::Rank2 &h_UU)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    chris_t out{};

    FOR (i, j, k)
    {
        out.LLL(i, j, k) = 0.5_rt * (d1_metric(j, i, k) + d1_metric(k, i, j) -
                                  d1_metric(j, k, i));
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

} // namespace TensorAlgebra

#endif /* TENSORALGEBRA_HPP_ */
