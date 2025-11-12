/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef TENSORALGEBRA_HPP_
#define TENSORALGEBRA_HPP_

// GRTeclyn includes
#include "AlwaysInline.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"

// AMReX includes
#include <AMReX_Array.H>
#include <AMReX_REAL.H>

// System includes
#include <array>

// This is a convenience function used to calculate the index that goes into
// Levi-Civita symbol which is actually represented as an amrex::Array1D

// NOLINTBEGIN(bugprone-easily-swappable-parameters,
// readability-identifier-length)
constexpr int index4D(int i, int j, int k, int l)
{
    return 64 * i + 16 * j + 4 * k + l;
}
// NOLINTEND(bugprone-easily-swappable-parameters,
// readability-identifier-length)

struct chris_t
{
    Tensor<3, amrex::Real> ULL;        //!< standard christoffel symbols
    Tensor<3, amrex::Real> LLL;        //!< 3 lower indices
    Tensor<1, amrex::Real> contracted; //!< contracted christoffel
};


struct chris_array_t
{
    amrex::Array3D<amrex::Real, 0, 3, 0, 3, 0, 3>
        ULL; //!< standard christoffel symbols
    amrex::Array3D<amrex::Real, 0, 3, 0, 3, 0, 3> LLL; //!< 3 lower indices
    amrex::Array1D<amrex::Real, 0, 3> contracted; //!< contracted christoffel
};

namespace TensorAlgebra
{
//  using Tensor<2, amrex::Real, 3> = amrex::Array2D<amrex::Real, 0, 3, 0, 3>

/// Computes determinant of a symmetric 3x3 matrix
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_determinant_sym(const amrex::Array2D<amrex::Real, 0, 3, 0, 3> &matrix)
{
    amrex::Real det = matrix(0, 0) * matrix(1, 1) * matrix(2, 2) +
                      2 * matrix(0, 1) * matrix(0, 2) * matrix(1, 2) -
                      matrix(0, 0) * matrix(1, 2) * matrix(1, 2) -
                      matrix(1, 1) * matrix(0, 2) * matrix(0, 2) -
                      matrix(2, 2) * matrix(0, 1) * matrix(0, 1);

    return det;
}

/// Computes determinant of a symmetric 3x3 matrix
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_determinant_sym(const Tensor<2, amrex::Real, 3> &matrix)
{
    amrex::Real det = matrix[0][0] * matrix[1][1] * matrix[2][2] +
                      2 * matrix[0][1] * matrix[0][2] * matrix[1][2] -
                      matrix[0][0] * matrix[1][2] * matrix[1][2] -
                      matrix[1][1] * matrix[0][2] * matrix[0][2] -
                      matrix[2][2] * matrix[0][1] * matrix[0][1];

    return det;
}

/// Computes the determinant of a general 3x3 matrix.
/// Note: for a symmetric matrix use the simplified function
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_determinant(const amrex::Array2D<amrex::Real, 0, 3, 0, 3> &matrix)
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

[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_determinant(const Tensor<2, amrex::Real, 3> &matrix)
{
    amrex::Real det =
        matrix[0][0] *
            (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
        matrix[0][1] *
            (matrix[2][2] * matrix[1][0] - matrix[1][2] * matrix[2][0]) +
        matrix[0][2] *
            (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
    return det;
}

/// Computes the inverse of a symmetric 3x3 matrix directly.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor<2, amrex::Real>
compute_inverse_sym(const Tensor<2, amrex::Real, 3> &matrix)
{
    amrex::Real deth         = compute_determinant_sym(matrix);
    amrex::Real deth_inverse = 1. / deth;
    Tensor<2, amrex::Real> h_UU;
    h_UU[0][0] = (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[1][2]) *
                 deth_inverse;
    h_UU[0][1] = (matrix[0][2] * matrix[1][2] - matrix[0][1] * matrix[2][2]) *
                 deth_inverse;
    h_UU[0][2] = (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]) *
                 deth_inverse;
    h_UU[1][1] = (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[0][2]) *
                 deth_inverse;
    h_UU[1][2] = (matrix[0][1] * matrix[0][2] - matrix[0][0] * matrix[1][2]) *
                 deth_inverse;
    h_UU[2][2] = (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[0][1]) *
                 deth_inverse;
    h_UU[1][0] = h_UU[0][1];
    h_UU[2][0] = h_UU[0][2];
    h_UU[2][1] = h_UU[1][2];

    return h_UU;
}

/// Computes the inverse of a symmetric 3x3 matrix directly.
[[nodiscard]] AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE amrex::Array2D<amrex::Real, 0, 3, 0, 3>
    compute_inverse_sym(const amrex::Array2D<amrex::Real, 0, 3, 0, 3> &matrix)
{
    amrex::Real deth         = compute_determinant_sym(matrix);
    amrex::Real deth_inverse = 1. / deth;
    amrex::Array2D<amrex::Real, 0, 3, 0, 3> h_UU{};
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
    h_UU(1, 0) = h_UU(0, 1);
    h_UU(2, 0) = h_UU(0, 2);
    h_UU(2, 1) = h_UU(1, 2);

    return h_UU;
}

/// Computes the inverse of a general 3x3 matrix.
/// Note: for a symmetric matrix use the simplified function
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor<2, amrex::Real>
compute_inverse(const Tensor<2, amrex::Real, 3> &matrix)
{
    amrex::Real deth         = compute_determinant(matrix);
    amrex::Real deth_inverse = 1. / deth;
    Tensor<2, amrex::Real> h_UU{};
    h_UU[0][0] = (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) *
                 deth_inverse;
    h_UU[1][1] = (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]) *
                 deth_inverse;
    h_UU[2][2] = (matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1]) *
                 deth_inverse;
    h_UU[1][0] = (matrix[2][0] * matrix[1][2] - matrix[1][0] * matrix[2][2]) *
                 deth_inverse;
    h_UU[0][1] = (matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2]) *
                 deth_inverse;
    h_UU[2][0] = (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]) *
                 deth_inverse;
    h_UU[0][2] = (matrix[0][1] * matrix[1][2] - matrix[1][1] * matrix[0][2]) *
                 deth_inverse;
    h_UU[2][1] = (matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1]) *
                 deth_inverse;
    h_UU[1][2] = (matrix[1][0] * matrix[0][2] - matrix[0][0] * matrix[1][2]) *
                 deth_inverse;

    return h_UU;
}

/// Computes the inverse of a general 3x3 matrix.
/// Note: for a symmetric matrix use the simplified function
[[nodiscard]] AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE amrex::Array2D<amrex::Real, 0, 3, 0, 3>
    compute_inverse(const amrex::Array2D<amrex::Real, 0, 3, 0, 3> &matrix)
{
    amrex::Real deth         = compute_determinant(matrix);
    amrex::Real deth_inverse = 1. / deth;
    amrex::Array2D<amrex::Real, 0, 3, 0, 3> h_UU{};
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
compute_trace(const amrex::Array2D<amrex::Real, 0, 3, 0, 3> &tensor_LL,
              const amrex::Array2D<amrex::Real, 0, 3, 0, 3> &inverse_metric)
{
    amrex::Real trace = 0.;
    FOR (i, j)
    {
        trace += inverse_metric(i, j) * tensor_LL(i, j);
    }
    return trace;
}

/// Computes the trace of a 2-Tensor with lower indices given an inverse metric.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_trace(const Tensor<2, amrex::Real> &tensor_LL,
              const Tensor<2, amrex::Real> &inverse_metric)
{
    amrex::Real trace = 0.;
    FOR (i, j)
    {
        trace += inverse_metric[i][j] * tensor_LL[i][j];
    }
    return trace;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_trace(const Tensor<2, amrex::Real> &tensor_LL,
              const amrex::Array1D<amrex::Real, 0, 6> &inverse_metric_sym)
{
    amrex::Real trace = 0.;
    FOR (i, j)
    {
        int idx  = i + j + ((i * j != 0) ? 1 : 0);
        trace   += inverse_metric_sym(idx) * tensor_LL[i][j];
    }
    return trace;
}

/// Computes the trace of a 1,1 Tensor (a matrix) - no metric required.
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_trace(const amrex::Array2D<amrex::Real, 0, 3, 0, 3> &tensor_UL)
{
    amrex::Real trace = 0.;
    FOR (i)
        trace += tensor_UL(i, i);
    return trace;
}

/// Computes the trace of a 1,1 Tensor (a matrix) - no metric required.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_trace(const Tensor<2, amrex::Real> &tensor_UL)
{
    amrex::Real trace = 0.;
    FOR (i)
        trace += tensor_UL[i][i];
    return trace;
}

/// I think this is right... A 1D array of 1D arrays.
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_trace(
    const amrex::Array1D<amrex::Array1D<amrex::Real, 0, 3>, 0, 3> &tensor_UL)
{
    amrex::Real trace = 0.;
    FOR (i)
        trace += tensor_UL(i)(i);
    return trace;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_trace(const Tensor<1, Tensor<1, amrex::Real>> &tensor_UL)
{
    amrex::Real trace = 0.;
    FOR (i)
        trace += tensor_UL[i][i];
    return trace;
}

/// Computes dot product of a vector and a covector (no metric required)
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_dot_product(const amrex::Array1D<amrex::Real, 0, 3> &vector_U,
                    const amrex::Array1D<amrex::Real, 0, 3> &covector_L)
{
    amrex::Real dot_product = 0.;
    FOR (i)
        dot_product += vector_U(i) * covector_L(i);
    return dot_product;
}

/// Computes dot product of a vector and a covector (no metric required)
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_dot_product(const amrex::Array1D<amrex::Real, 0, 3> &vector_U,
                    const Tensor<1, amrex::Real> &covector_L)
{
    amrex::Real dot_product = 0.;
    FOR (i)
        dot_product += vector_U(i) * covector_L[i];
    return dot_product;
}

/// Computes dot product of a vector and a covector (no metric required)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_dot_product(const Tensor<1, amrex::Real> &vector_U,
                    const Tensor<1, amrex::Real> &covector_L)
{
    amrex::Real dot_product = 0.;
    FOR (i)
        dot_product += vector_U[i] * covector_L[i];
    return dot_product;
}

/// Computes dot product of two covectors given an inverse metric or
/// the dot product of two vectors given a metric.
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_dot_product(
    const amrex::Array1D<amrex::Real, 0, 3> &covector1_L,
    const amrex::Array1D<amrex::Real, 0, 3> &covector2_L,
    const amrex::Array2D<amrex::Real, 0, 3, 0, 3> &inverse_metric)
{
    amrex::Real dot_product = 0.;
    FOR (m, n)
    {
        dot_product += inverse_metric(m, n) * covector1_L(m) * covector2_L(n);
    }
    return dot_product;
}

/// Computes dot product of two covectors given an inverse metric or
/// the dot product of two vectors given a metric.
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_dot_product(
    const Tensor<1, amrex::Real> &covector1_L,
    const Tensor<1, amrex::Real> &covector2_L,
    const amrex::Array2D<amrex::Real, 0, 3, 0, 3> &inverse_metric)
{
    amrex::Real dot_product = 0.;
    FOR (m, n)
    {
        dot_product += inverse_metric(m, n) * covector1_L[m] * covector2_L[n];
    }
    return dot_product;
}

/// Computes dot product of two covectors given an inverse metric or
/// the dot product of two vectors given a metric.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_dot_product(const Tensor<1, amrex::Real> &covector1_L,
                    const Tensor<1, amrex::Real> &covector2_L,
                    const Tensor<2, amrex::Real> &inverse_metric)
{
    amrex::Real dot_product = 0.;
    FOR (m, n)
    {
        dot_product += inverse_metric[m][n] * covector1_L[m] * covector2_L[n];
    }
    return dot_product;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_dot_product(const Tensor<1, amrex::Real> &covector1_L,
                    const Tensor<1, amrex::Real> &covector2_L,
                    const amrex::Array1D<amrex::Real, 0, 6> &inverse_metric_sym)
{
    amrex::Real dot_product = 0.;
    FOR (m, n)
    {
        int idx = m + n + ((m * n != 0) ? 1 : 0);
        dot_product +=
            inverse_metric_sym(idx) * covector1_L[m] * covector2_L[n];
    }
    return dot_product;
}

/// Removes the trace of a 2-Tensor with lower indices given a metric and an
/// inverse metric.  Or a Tensor with upper indices given an inverse metric and
/// a metric.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
// NOLINTBEGIN(bugprone-easily-swappable-parameters)
make_trace_free(amrex::Array2D<amrex::Real, 0, 3, 0, 3> &tensor_LL,
                const amrex::Array2D<amrex::Real, 0, 3, 0, 3> &metric,
                const amrex::Array2D<amrex::Real, 0, 3, 0, 3> &inverse_metric)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    auto trace                  = compute_trace(tensor_LL, inverse_metric);
    double one_over_gr_spacedim = 1. / ((double)GR_SPACEDIM);
    FOR (i, j)
    {
        tensor_LL(i, j) -= one_over_gr_spacedim * metric(i, j) * trace;
    }
}

/// Removes the trace of a 2-Tensor with lower indices given a metric and an
/// inverse metric.  Or a Tensor with upper indices given an inverse metric and
/// a metric.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
// NOLINTBEGIN(bugprone-easily-swappable-parameters)
make_trace_free(Tensor<2, amrex::Real> &tensor_LL,
                const Tensor<2, amrex::Real> &metric,
                const Tensor<2, amrex::Real> &inverse_metric)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    auto trace                  = compute_trace(tensor_LL, inverse_metric);
    double one_over_gr_spacedim = 1. / ((double)GR_SPACEDIM);
    FOR (i, j)
    {
        tensor_LL[i][j] -= one_over_gr_spacedim * metric[i][j] * trace;
    }
}

/// Makes a 2-Tensor symmetric
template <int size>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
make_symmetric(amrex::Array2D<amrex::Real, 0, size, 0, size> &tensor_LL)
{
    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < i; ++j)
        {
            tensor_LL(i, j) = 0.5 * (tensor_LL(i, j) + tensor_LL(j, i));
            tensor_LL(j, i) = tensor_LL(i, j);
        }
    }
}

/// Makes a 2-Tensor symmetric
template <int size>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
make_symmetric(Tensor<2, amrex::Real, size> &tensor_LL)
{
    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < i; ++j)
        {
            tensor_LL[i][j] = 0.5 * (tensor_LL[i][j] + tensor_LL[j][i]);
            tensor_LL[j][i] = tensor_LL[i][j];
        }
    }
}

/// Raises the index of a covector
[[nodiscard]] AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE amrex::Array1D<amrex::Real, 0, 3>
    raise_all(const amrex::Array1D<amrex::Real, 0, 3> &tensor_L,
              const amrex::Array2D<amrex::Real, 0, 3, 0, 3> &inverse_metric)
{
    amrex::Array1D<amrex::Real, 0, 3> tensor_U{};

    FOR (i)
    {
        tensor_U(i) = 0.;
    }

    FOR (i, j)
    {
        tensor_U(i) += inverse_metric(i, j) * tensor_L(j);
    }
    return tensor_U;
}

/// Raises the index of a covector
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor<1, amrex::Real>
raise_all(const Tensor<1, amrex::Real> &tensor_L,
          const Tensor<2, amrex::Real> &inverse_metric)
{
    Tensor<1, amrex::Real> tensor_U = 0.;
    FOR (i, j)
    {
        tensor_U[i] += inverse_metric[i][j] * tensor_L[j];
    }
    return tensor_U;
}

/// Raises the indices of a 2-Tensor
[[nodiscard]] AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE amrex::Array2D<amrex::Real, 0, 3, 0, 3>
    raise_all(const amrex::Array2D<amrex::Real, 0, 3, 0, 3> &tensor_LL,
              const amrex::Array2D<amrex::Real, 0, 3, 0, 3> &inverse_metric)
{
    amrex::Array2D<amrex::Real, 0, 3, 0, 3> tensor_UU{};

    FOR (i, j)
    {
        tensor_UU(i, j) = 0.;
    }

    FOR (i, j, k, l)
    {
        tensor_UU(i, j) +=
            inverse_metric(i, k) * inverse_metric(j, l) * tensor_LL(k, l);
    }
    return tensor_UU;
}

/// Raises the indices of a 2-Tensor
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor<2, amrex::Real>
raise_all(const Tensor<2, amrex::Real> &tensor_LL,
          const Tensor<2, amrex::Real> &inverse_metric)
{
    Tensor<2, amrex::Real> tensor_UU = 0.;
    FOR (i, j, k, l)
    {
        tensor_UU[i][j] +=
            inverse_metric[i][k] * inverse_metric[j][l] * tensor_LL[k][l];
    }
    return tensor_UU;
}

/// Lowers the indices of a vector
/// Note: same functionality as raise; included to improve readability
[[nodiscard]] AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE amrex::Array1D<amrex::Real, 0, 3>
    lower_all(const amrex::Array1D<amrex::Real, 0, 3> &tensor_U,
              const amrex::Array2D<amrex::Real, 0, 3, 0, 3> &metric)
{ // The code for lowering is exactly the same as for raising
    return raise_all(tensor_U, metric);
}

/// Lowers the indices of a 2-Tensor
/// Note: same functionality as raise; included to improve readability
[[nodiscard]] AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE amrex::Array2D<amrex::Real, 0, 3, 0, 3>
    lower_all(const amrex::Array2D<amrex::Real, 0, 3, 0, 3> &tensor_UU,
              const amrex::Array2D<amrex::Real, 0, 3, 0, 3> &metric)
{ // The code for lowering is exactly the same as for raising
    return raise_all(tensor_UU, metric);
}

/// Lowers the indices of a vector
/// Note: same functionality as raise; included to improve readibility
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor<1, amrex::Real>
lower_all(const Tensor<1, amrex::Real> &tensor_U,
          const Tensor<2, amrex::Real> &metric)
{ // The code for lowering is exactly the same as for raising
    return raise_all(tensor_U, metric);
}

/// Lowers the indices of a 2-Tensor
/// Note: same functionality as raise; included to improve readibility
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor<2, amrex::Real>
lower_all(const Tensor<2, amrex::Real> &tensor_UU,
          const Tensor<2, amrex::Real> &metric)
{ // The code for lowering is exactly the same as for raising
    return raise_all(tensor_UU, metric);
}

/// Computes the (i,j) component of the Kronecker delta
constexpr int delta(int i, int j) { return static_cast<int>(i == j); }

/// Computes the levi-civita symbol (3D, NB, symbol, not the Tensor)
[[nodiscard]] AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE amrex::Array3D<amrex::Real, 0, 3, 0, 3, 0, 3>
    epsilon_array()
{
    amrex::Array3D<amrex::Real, 0, 3, 0, 3, 0, 3> epsilon{};

    FOR (i, j, k)
    {
        epsilon(i, j, k) = 0.;
    }

    epsilon(0, 1, 2) = 1.0;
    epsilon(1, 2, 0) = 1.0;
    epsilon(2, 0, 1) = 1.0;
    epsilon(0, 2, 1) = -1.0;
    epsilon(2, 1, 0) = -1.0;
    epsilon(1, 0, 2) = -1.0;

    return epsilon;
}

/// Computes the levi-civita symbol (3D, NB, symbol, not the Tensor)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor<3, double> epsilon()
{
    Tensor<3, double> epsilon = {0.};
    epsilon[0][1][2]          = 1.0;
    epsilon[1][2][0]          = 1.0;
    epsilon[2][0][1]          = 1.0;
    epsilon[0][2][1]          = -1.0;
    epsilon[2][1][0]          = -1.0;
    epsilon[1][0][2]          = -1.0;

    return epsilon;
}

/// Computes the levi-civita symbol (4D, NB, symbol, not the Tensor)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor<4, double, 4> epsilon4D()
{
    Tensor<4, double, 4> epsilon4D = {0.0};
    epsilon4D[0][1][2][3]          = 1.0;
    epsilon4D[0][1][3][2]          = -1.0;
    epsilon4D[0][3][1][2]          = 1.0;
    epsilon4D[0][3][2][1]          = -1.0;
    epsilon4D[0][2][1][3]          = -1.0;
    epsilon4D[0][2][3][1]          = 1.0;

    epsilon4D[1][0][2][3] = -1.0;
    epsilon4D[1][2][0][3] = 1.0;
    epsilon4D[1][2][3][0] = -1.0;
    epsilon4D[1][3][2][0] = 1.0;
    epsilon4D[1][3][0][2] = -1.0;
    epsilon4D[1][0][3][2] = 1.0;

    epsilon4D[2][0][1][3] = 1.0;
    epsilon4D[2][0][3][1] = -1.0;
    epsilon4D[2][3][0][1] = 1.0;
    epsilon4D[2][3][1][0] = -1.0;
    epsilon4D[2][1][3][0] = 1.0;
    epsilon4D[2][1][0][3] = -1.0;

    epsilon4D[3][0][1][2] = -1.0;
    epsilon4D[3][1][0][2] = 1.0;
    epsilon4D[3][1][2][0] = -1.0;
    epsilon4D[3][2][1][0] = 1.0;
    epsilon4D[3][2][0][1] = -1.0;
    epsilon4D[3][0][2][1] = 1.0;

    return epsilon4D;
}

/// Computes the levi-civita symbol (4D, NB, symbol, not the Tensor)
[[nodiscard]] AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE amrex::Array1D<amrex::Real, 0, 256>
    epsilon4D_array()
{
    amrex::Array1D<amrex::Real, 0, 256> epsilon4D{};

    for (unsigned int i = 0; i < 256; i++)
    {
        epsilon4D(i) = 0.0;
    }

    // Fortran order!
    epsilon4D(index4D(0, 1, 2, 3)) = 1.0;
    epsilon4D(index4D(0, 1, 3, 2)) = -1.0;
    epsilon4D(index4D(0, 3, 1, 2)) = 1.0;
    epsilon4D(index4D(0, 3, 2, 1)) = -1.0;
    epsilon4D(index4D(0, 2, 1, 3)) = -1.0;
    epsilon4D(index4D(0, 2, 3, 1)) = 1.0;

    epsilon4D(index4D(1, 0, 2, 3)) = -1.0;
    epsilon4D(index4D(1, 2, 0, 3)) = 1.0;
    epsilon4D(index4D(1, 2, 3, 0)) = -1.0;
    epsilon4D(index4D(1, 3, 2, 0)) = 1.0;
    epsilon4D(index4D(1, 3, 0, 2)) = -1.0;
    epsilon4D(index4D(1, 0, 3, 2)) = 1.0;

    epsilon4D(index4D(2, 0, 1, 3)) = 1.0;
    epsilon4D(index4D(2, 0, 3, 1)) = -1.0;
    epsilon4D(index4D(2, 3, 0, 1)) = 1.0;
    epsilon4D(index4D(2, 3, 1, 0)) = -1.0;
    epsilon4D(index4D(2, 1, 3, 0)) = 1.0;
    epsilon4D(index4D(2, 1, 0, 3)) = -1.0;

    epsilon4D(index4D(3, 0, 1, 2)) = -1.0;
    epsilon4D(index4D(3, 1, 0, 2)) = 1.0;
    epsilon4D(index4D(3, 1, 2, 0)) = -1.0;
    epsilon4D(index4D(3, 2, 1, 0)) = 1.0;
    epsilon4D(index4D(3, 2, 0, 1)) = -1.0;
    epsilon4D(index4D(3, 0, 2, 1)) = 1.0;

    return epsilon4D;
}

/// Computes the conformal christoffel symbol
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE chris_t
compute_christoffel(const Tensor<2, Tensor<1, amrex::Real>> &d1_metric,
                    const Tensor<2, amrex::Real> &h_UU)
{
    chris_t out{};

    FOR (i, j, k)
    {
        out.LLL[i][j][k] = 0.5 * (d1_metric[j][i][k] + d1_metric[k][i][j] -
                                  d1_metric[j][k][i]);
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

/// Computes the conformal christoffel symbol
[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE chris_array_t
compute_christoffel(const amrex::Array2D<amrex::Array1D<amrex::Real, 0, 3>, 0,
                                         3, 0, 3> &d1_metric,
                    const amrex::Array2D<amrex::Real, 0, 3, 0, 3> &h_UU)
{
    chris_array_t out{};

    FOR (i, j, k)
    {
        out.LLL(i, j, k) = 0.5 * (d1_metric(j, i)(k) + d1_metric(k, i)(j) -
                                  d1_metric(j, k)(i));
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
