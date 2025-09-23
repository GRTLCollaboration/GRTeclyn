/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef TENSORALGEBRA_HPP_
#define TENSORALGEBRA_HPP_

#include "AlwaysInline.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include <AMReX_REAL.H>
#include <array>

struct chris_t
{
    Tensor<3, amrex::Real> ULL;        //!< standard christoffel symbols
    Tensor<3, amrex::Real> LLL;        //!< 3 lower indices
    Tensor<1, amrex::Real> contracted; //!< contracted christoffel
};

namespace TensorAlgebra
{
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
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
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

/// Computes the inverse of a general 3x3 matrix.
/// Note: for a symmetric matrix use the simplified function
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor<2, amrex::Real>
compute_inverse(const Tensor<2, amrex::Real, 3> &matrix)
{
    amrex::Real deth         = compute_determinant(matrix);
    amrex::Real deth_inverse = 1. / deth;
    Tensor<2, amrex::Real> h_UU;
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

/// Computes the trace of a 1,1 Tensor (a matrix) - no metric required.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_trace(const Tensor<2, amrex::Real> &tensor_UL)
{
    amrex::Real trace = 0.;
    FOR (i)
        trace += tensor_UL[i][i];
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
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE chris_t compute_christoffel(
    const Tensor<3, amrex::Real> &d1_metric, const Tensor<2, amrex::Real> &h_UU)
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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor<3, amrex::Real>
compute_phys_chris(const Tensor<1, amrex::Real> &d1_chi,
                   const amrex::Real &vars_chi,
                   const Tensor<2, amrex::Real> &vars_h,
                   const Tensor<2, amrex::Real> &h_UU,
                   const Tensor<3, amrex::Real> &chris_ULL)
{
    Tensor<3, amrex::Real> chris_phys;
    FOR (i, j, k)
    {
        chris_phys[i][j][k] =
            chris_ULL[i][j][k] -
            0.5 / vars_chi *
                (delta(i, k) * d1_chi[j] + delta(i, j) * d1_chi[k]);
        FOR (m)
        {
            chris_phys[i][j][k] +=
                0.5 / vars_chi * vars_h[j][k] * h_UU[i][m] * d1_chi[m];
        }
    }
    return chris_phys;
}
} // namespace TensorAlgebra

#endif /* TENSORALGEBRA_HPP_ */
