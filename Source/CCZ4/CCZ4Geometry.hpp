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
#include "FourthOrderDerivatives.hpp"
#include "TensorAlgebra.hpp"

//! A structure for the decomposed elements of the Energy Momentum Tensor in
//! 3+1D
struct emtensor_t
{
    amrex::Array2D<amrex::Real, 0, 3, 0, 3> S; //!< S_ij = T_ij
    amrex::Array1D<amrex::Real, 0, 3> j;       //!< j_i = T_ia_n^a
    amrex::Real trS;                           //!< trS = S^i_i
    amrex::Real rho;                           //!< rho = T_ab n^a n^b
};

struct ricci_t
{
    amrex::Array2D<amrex::Real, 0, 3, 0, 3> LL; // Ricci with two indices down
    amrex::Real scalar{};                       // Ricci scalar
};

struct ricci_t_array
{
    amrex::Array2D<amrex::Real, 0, 3, 0, 3> LL; // Ricci with two indices down
    amrex::Real scalar{};                       // Ricci scalar
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

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_metric_determinant(const amrex::CellData<const amrex::Real> &h)
{
    amrex::Real det = h[var_idx(c_h11, 0, 0)] * h[var_idx(c_h11, 1, 1)] *
                          h[var_idx(c_h11, 2, 2)] +
                      2 * h[var_idx(c_h11, 0, 1)] * h[var_idx(c_h11, 0, 2)] *
                          h[var_idx(c_h11, 1, 2)] -
                      h[var_idx(c_h11, 0, 0)] * h[var_idx(c_h11, 1, 2)] *
                          h[var_idx(c_h11, 1, 2)] -
                      h[var_idx(c_h11, 1, 1)] * h[var_idx(c_h11, 0, 2)] *
                          h[var_idx(c_h11, 0, 2)] -
                      h[var_idx(c_h11, 2, 2)] * h[var_idx(c_h11, 0, 1)] *
                          h[var_idx(c_h11, 0, 1)];

    return det;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Array2D<amrex::Real, 0, 3, 0, 3>
compute_inverse_metric(const CCZ4Vars &vars)
{
    amrex::Real det_h         = compute_metric_determinant(vars);
    amrex::Real det_h_inverse = 1. / det_h;
    amrex::Array2D<amrex::Real, 0, 3, 0, 3> h_UU{};

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

AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor<2, amrex::Real>
compute_inverse_metric(const amrex::CellData<const amrex::Real> &h)
{
    amrex::Real det_h         = compute_metric_determinant(h);
    amrex::Real det_h_inverse = 1. / det_h;
    Tensor<2, amrex::Real> h_UU;

    h_UU[0][0] = (h[var_idx(c_h11, 1, 1)] * h[var_idx(c_h11, 2, 2)] -
                  h[var_idx(c_h11, 1, 2)] * h[var_idx(c_h11, 2, 1)]) *
                 det_h_inverse;
    h_UU[0][1] = (h[var_idx(c_h11, 2, 0)] * h[var_idx(c_h11, 1, 2)] -
                  h[var_idx(c_h11, 1, 0)] * h[var_idx(c_h11, 2, 2)]) *
                 det_h_inverse;
    h_UU[0][2] = (h[var_idx(c_h11, 1, 0)] * h[var_idx(c_h11, 2, 1)] -
                  h[var_idx(c_h11, 2, 0)] * h[var_idx(c_h11, 1, 1)]) *
                 det_h_inverse;
    h_UU[1][1] = (h[var_idx(c_h11, 0, 0)] * h[var_idx(c_h11, 2, 2)] -
                  h[var_idx(c_h11, 2, 0)] * h[var_idx(c_h11, 0, 2)]) *
                 det_h_inverse;
    h_UU[1][2] = (h[var_idx(c_h11, 0, 1)] * h[var_idx(c_h11, 2, 0)] -
                  h[var_idx(c_h11, 0, 0)] * h[var_idx(c_h11, 2, 1)]) *
                 det_h_inverse;
    h_UU[2][2] = (h[var_idx(c_h11, 0, 0)] * h[var_idx(c_h11, 1, 1)] -
                  h[var_idx(c_h11, 0, 1)] * h[var_idx(c_h11, 1, 0)]) *
                 det_h_inverse;
    h_UU[1][0] = h_UU[0][1];
    h_UU[2][0] = h_UU[0][2];
    h_UU[2][1] = h_UU[1][2];

    return h_UU;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor<1, amrex::Real, 6>
compute_inverse_metric_sym_tensor(const amrex::CellData<const amrex::Real> &h)
{
    amrex::Real det_h         = compute_metric_determinant(h);
    amrex::Real det_h_inverse = 1. / det_h;
    Tensor<1, amrex::Real, 6> h_UU;
    h_UU[0] = (h[var_idx(c_h11, 1, 1)] * h[var_idx(c_h11, 2, 2)] -
               h[var_idx(c_h11, 1, 2)] * h[var_idx(c_h11, 2, 1)]) *
              det_h_inverse;
    h_UU[1] = (h[var_idx(c_h11, 2, 0)] * h[var_idx(c_h11, 1, 2)] -
               h[var_idx(c_h11, 1, 0)] * h[var_idx(c_h11, 2, 2)]) *
              det_h_inverse;
    h_UU[2] = (h[var_idx(c_h11, 1, 0)] * h[var_idx(c_h11, 2, 1)] -
               h[var_idx(c_h11, 2, 0)] * h[var_idx(c_h11, 1, 1)]) *
              det_h_inverse;
    h_UU[3] = (h[var_idx(c_h11, 0, 0)] * h[var_idx(c_h11, 2, 2)] -
               h[var_idx(c_h11, 2, 0)] * h[var_idx(c_h11, 0, 2)]) *
              det_h_inverse;
    h_UU[4] = (h[var_idx(c_h11, 0, 1)] * h[var_idx(c_h11, 2, 0)] -
               h[var_idx(c_h11, 0, 0)] * h[var_idx(c_h11, 2, 1)]) *
              det_h_inverse;
    h_UU[5] = (h[var_idx(c_h11, 0, 0)] * h[var_idx(c_h11, 1, 1)] -
               h[var_idx(c_h11, 0, 1)] * h[var_idx(c_h11, 1, 0)]) *
              det_h_inverse;

    return h_UU;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Array1D<amrex::Real, 0, 6>
compute_inverse_metric_sym(const amrex::CellData<const amrex::Real> &h)
{
    amrex::Real det_h         = compute_metric_determinant(h);
    amrex::Real det_h_inverse = 1. / det_h;
    amrex::Array1D<amrex::Real, 0, 6> h_UU;
    h_UU(0) = (h[var_idx(c_h11, 1, 1)] * h[var_idx(c_h11, 2, 2)] -
               h[var_idx(c_h11, 1, 2)] * h[var_idx(c_h11, 2, 1)]) *
              det_h_inverse;
    h_UU(1) = (h[var_idx(c_h11, 2, 0)] * h[var_idx(c_h11, 1, 2)] -
               h[var_idx(c_h11, 1, 0)] * h[var_idx(c_h11, 2, 2)]) *
              det_h_inverse;
    h_UU(2) = (h[var_idx(c_h11, 1, 0)] * h[var_idx(c_h11, 2, 1)] -
               h[var_idx(c_h11, 2, 0)] * h[var_idx(c_h11, 1, 1)]) *
              det_h_inverse;
    h_UU(3) = (h[var_idx(c_h11, 0, 0)] * h[var_idx(c_h11, 2, 2)] -
               h[var_idx(c_h11, 2, 0)] * h[var_idx(c_h11, 0, 2)]) *
              det_h_inverse;
    h_UU(4) = (h[var_idx(c_h11, 0, 1)] * h[var_idx(c_h11, 2, 0)] -
               h[var_idx(c_h11, 0, 0)] * h[var_idx(c_h11, 2, 1)]) *
              det_h_inverse;
    h_UU(5) = (h[var_idx(c_h11, 0, 0)] * h[var_idx(c_h11, 1, 1)] -
               h[var_idx(c_h11, 0, 1)] * h[var_idx(c_h11, 1, 0)]) *
              det_h_inverse;

    return h_UU;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Array2D<amrex::Real, 0, 3, 0, 3>
compute_inverse_metric_array(const amrex::CellData<const amrex::Real> &h)
{
    amrex::Real det_h         = compute_metric_determinant(h);
    amrex::Real det_h_inverse = 1. / det_h;
    amrex::Array2D<amrex::Real, 0, 3, 0, 3> h_UU;
    h_UU(0, 0) = (h[var_idx(c_h11, 1, 1)] * h[var_idx(c_h11, 2, 2)] -
                  h[var_idx(c_h11, 1, 2)] * h[var_idx(c_h11, 2, 1)]) *
                 det_h_inverse;
    h_UU(0, 1) = (h[var_idx(c_h11, 2, 0)] * h[var_idx(c_h11, 1, 2)] -
                  h[var_idx(c_h11, 1, 0)] * h[var_idx(c_h11, 2, 2)]) *
                 det_h_inverse;
    h_UU(0, 2) = (h[var_idx(c_h11, 1, 0)] * h[var_idx(c_h11, 2, 1)] -
                  h[var_idx(c_h11, 2, 0)] * h[var_idx(c_h11, 1, 1)]) *
                 det_h_inverse;
    h_UU(1, 1) = (h[var_idx(c_h11, 0, 0)] * h[var_idx(c_h11, 2, 2)] -
                  h[var_idx(c_h11, 2, 0)] * h[var_idx(c_h11, 0, 2)]) *
                 det_h_inverse;
    h_UU(1, 2) = (h[var_idx(c_h11, 0, 1)] * h[var_idx(c_h11, 2, 0)] -
                  h[var_idx(c_h11, 0, 0)] * h[var_idx(c_h11, 2, 1)]) *
                 det_h_inverse;
    h_UU(2, 2) = (h[var_idx(c_h11, 0, 0)] * h[var_idx(c_h11, 1, 1)] -
                  h[var_idx(c_h11, 0, 1)] * h[var_idx(c_h11, 1, 0)]) *
                 det_h_inverse;
    h_UU(1, 0) = h_UU(0, 1);
    h_UU(2, 0) = h_UU(0, 2);
    h_UU(2, 1) = h_UU(1, 2);

    return h_UU;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_trace_A(const CCZ4Vars &vars)
{
    amrex::Array2D<amrex::Real, 0, 3, 0, 3> inverse_metric =
        compute_inverse_metric(vars);
    amrex::Real trace_A = 0.0;
    FOR (i, j)
    {
        trace_A += inverse_metric(i, j) * vars.A(i, j);
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

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Array2D<amrex::Real, 0, 3, 0, 3>
compute_A_UU(const CCZ4Vars &vars,
             const amrex::Array2D<amrex::Real, 0, 3, 0, 3> &inverse_metric)
{
    amrex::Array2D<amrex::Real, 0, 3, 0, 3> A_UU{};
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

AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor<2, amrex::Real>
compute_A_UU(const amrex::CellData<amrex::Real const> &cell_data,
             const Tensor<2, amrex::Real> &inverse_metric)
{
    Tensor<2, amrex::Real> A_UU;
    FOR (i, j)
    {
        A_UU[i][j] = 0.0;
        FOR (k, l)
        {
            A_UU[i][j] += inverse_metric[i][k] * inverse_metric[j][l] *
                          cell_data[var_idx(c_A11, k, l)];
        }
    }
    return A_UU;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Array2D<amrex::Real, 0, 3, 0, 3>
compute_A_UU(const amrex::CellData<amrex::Real const> &cell_data,
             const amrex::Array1D<amrex::Real, 0, 6> &inverse_metric_sym)
{
    amrex::Array2D<amrex::Real, 0, 3, 0, 3> A_UU;
    FOR (i, j)
    {
        A_UU(i, j) = 0.0;
        FOR (k, l)
        {
            A_UU(i, j) += inverse_metric_sym(SYMM_IDX(i, k)) *
                          inverse_metric_sym(SYMM_IDX(j, l)) *
                          cell_data[var_idx(c_A11, k, l)];
        }
    }
    return A_UU;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Array2D<amrex::Real, 0, 3, 0, 3>
compute_A_UU(const amrex::CellData<amrex::Real const> &cell_data,
             const amrex::Array2D<amrex::Real, 0, 3, 0, 3> &inverse_metric)
{
    amrex::Array2D<amrex::Real, 0, 3, 0, 3> A_UU;
    FOR (i, j)
    {
        A_UU(i, j) = 0.0;
        FOR (k, l)
        {
            A_UU(i, j) += inverse_metric(i, k) * inverse_metric(j, l) *
                          cell_data[var_idx(c_A11, k, l)];
        }
    }
    return A_UU;
}

// This is A_ij A^ij
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_Aij_squared(const CCZ4Vars &vars)
{
    amrex::Array2D<amrex::Real, 0, 3, 0, 3> inverse_metric =
        compute_inverse_metric(vars);
    amrex::Real Aij_squared = 0.0;
    FOR (i, j, k, l)
    {
        Aij_squared += inverse_metric(i, k) * inverse_metric(j, l) *
                       vars.A(i, j) * vars.A(k, l);
    }
    return Aij_squared;
}

// This is A_ij A^ij
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real compute_Aij_squared(
    const CCZ4Vars &vars,
    const amrex::Array2D<amrex::Real, 0, 3, 0, 3> &inverse_metric)
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
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real compute_Aij_squared_with_A_UU(
    const CCZ4Vars &vars, const amrex::Array2D<amrex::Real, 0, 3, 0, 3> &A_UU)
{
    amrex::Real Aij_squared = 0.0;
    FOR (i, j)
    {
        Aij_squared += A_UU(i, j) * vars.A(i, j);
    }
    return Aij_squared;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_Aij_squared(const amrex::CellData<const amrex::Real> &A,
                    const Tensor<2, amrex::Real> &inverse_metric)
{
    amrex::Real Aij_squared = 0.0;
    FOR (i, j, k, l)
    {
        Aij_squared += inverse_metric[i][k] * inverse_metric[j][l] *
                       A[var_idx(c_A11, i, j)] * A[var_idx(c_A11, k, l)];
    }
    return Aij_squared;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_Aij_squared(const amrex::CellData<const amrex::Real> &A,
                    const amrex::Array1D<amrex::Real, 0, 6> &inverse_metric_sym)
{
    amrex::Real Aij_squared = 0.0;
    FOR (i, j, k, l)
    {

        Aij_squared += inverse_metric_sym(SYMM_IDX(i, k)) *
                       inverse_metric_sym(SYMM_IDX(j, l)) *
                       A[var_idx(c_A11, i, j)] * A[var_idx(c_A11, k, l)];
    }
    return Aij_squared;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_Aij_squared(const amrex::CellData<const amrex::Real> &A,
                    const amrex::Array2D<amrex::Real, 0, 3, 0, 3> &A_UU)
{
    amrex::Real Aij_squared = 0.0;
    FOR (i, j)
    {
        Aij_squared += A[var_idx(c_A11, i, j)] * A_UU(i, j);
    }
    return Aij_squared;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_Aij_squared(const amrex::CellData<const amrex::Real> &A,
                    const Tensor<2, amrex::Real> &A_UU)
{
    amrex::Real Aij_squared = 0.0;
    FOR (i, j)
    {
        Aij_squared += A[var_idx(c_A11, i, j)] * A_UU[i][j];
    }
    return Aij_squared;
}

/// Computes the conformal christoffel symbol
AMREX_GPU_DEVICE AMREX_FORCE_INLINE chris_array_t compute_christoffel(
    const CCZ4D1Vars &d1, const amrex::Array2D<amrex::Real, 0, 3, 0, 3> &h_UU)
{
    chris_array_t out{};

    FOR (i, j, k)
    {
        out.LLL(i, j, k) =
            0.5 * (d1.h(j, i, k) + d1.h(k, i, j) - d1.h(j, k, i));
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

/// Computes the conformal christoffel symbol - using AMReX Arrays

AMREX_GPU_DEVICE AMREX_FORCE_INLINE chris_t_array
compute_christoffel(const amrex::Array2D<amrex::Real, 0, 6, 0, 3> &d1_h,
                    const amrex::Array1D<amrex::Real, 0, 6> &h_UU)
{
    chris_t_array out{};

    FOR (i, j, k)
    {
        out.LLL(i, j, k) =
            0.5 * (d1_h(SYMM_IDX(j, i), k) + d1_h(SYMM_IDX(k, i), j) -
                   d1_h(SYMM_IDX(j, k), i));
    }

    // FOR (i)
    // {
    //     out.LLL(i, i, i) = 0.5 * d1_h(SYMM_IDX(i, i), i);
    // }

    // out.LLL(0, 0, 1) = 0.5 * d1_h(0, 1);
    // out.LLL(0, 0, 2) = 0.5 * d1_h(0, 2);

    // out.LLL(0, 1, 0) = 0.5 * d1_h(0, 1);
    // out.LLL(0, 2, 0) = 0.5 * d1_h(0, 2);

    // out.LLL(0, 1, 1) = d1_h(1, 1) - 0.5 * d1_h(4, 0);
    // out.LLL(0, 2, 2) = d1_h(2, 2) - 0.5 * d1_h(5, 0);

    // out.LLL(1, 0, 0) = d1_h(1, 0) - 0.5 * d1_h(0, 1);
    // out.LLL(1, 2, 2) = d1_h(4, 2) - 0.5 * d1_h(5, 1);

    // out.LLL(1, 0, 1) = 0.5 * d1_h(3, 0);
    // out.LLL(1, 1, 0) = 0.5 * d1_h(3, 0);

    // out.LLL(1, 1, 2) = 0.5 * d1_h(3, 2);
    // out.LLL(1, 2, 1) = 0.5 * d1_h(3, 2);

    // out.LLL(2, 0, 0) = d1_h(2, 0) - 0.5 * d1_h(0, 2);
    // out.LLL(2, 1, 1) = d1_h(4, 1) - 0.5 * d1_h(3, 2);

    // out.LLL(2, 1, 2) = 0.5 * d1_h(5, 1);
    // out.LLL(2, 2, 1) = 0.5 * d1_h(5, 1);

    // out.LLL(2, 2, 0) = 0.5 * d1_h(5, 0);
    // out.LLL(2, 0, 2) = 0.5 * d1_h(5, 0);

    // ////These have all different indices

    // out.LLL(0, 1, 2) = 0.5 * (d1_h(1, 2) + d1_h(2, 0) + d1_h(4, 0));
    // out.LLL(0, 2, 1) = out.LLL(0, 1, 2);

    // out.LLL(1, 0, 2) = 0.5 * (d1_h(1, 2) + d1_h(4, 0) - d1_h(2, 0));
    // out.LLL(1, 2, 0) = out.LLL(1, 0, 2);

    // out.LLL(2, 0, 1) = 0.5 * (d1_h(3, 1) + d1_h(4, 0) - d1_h(1, 2));
    // out.LLL(2, 1, 0) = out.LLL(2, 0, 1);

    FOR (i, j)
    {
        out.ULL(0, i, j) = h_UU(0) * out.LLL(0, i, j) +
                           h_UU(1) * out.LLL(1, i, j) +
                           h_UU(2) * out.LLL(2, i, j);
        out.ULL(1, i, j) = h_UU(1) * out.LLL(0, i, j) +
                           h_UU(3) * out.LLL(1, i, j) +
                           h_UU(4) * out.LLL(2, i, j);
        out.ULL(2, i, j) = h_UU(2) * out.LLL(0, i, j) +
                           h_UU(4) * out.LLL(1, i, j) +
                           h_UU(5) * out.LLL(2, i, j);
    }

    // FOR (i, j, k)
    // {
    //     out.ULL(i, j, k) = 0;
    //     FOR (l)
    //     {
    //         out.ULL(i, j, k) += h_UU(SYMM_IDX(i, l)) * out.LLL(l, j, k);
    //     }
    // }
    FOR (i)
    {
        out.contracted(i) =
            h_UU(0) * out.ULL(i, 0, 0) + h_UU(1) * out.ULL(i, 0, 1) +
            h_UU(2) * out.ULL(i, 0, 2) + h_UU(1) * out.ULL(i, 1, 0) +
            h_UU(3) * out.ULL(i, 1, 1) + h_UU(4) * out.ULL(i, 1, 2) +
            h_UU(2) * out.ULL(i, 2, 0) + h_UU(4) * out.ULL(i, 2, 1) +
            h_UU(5) * out.ULL(i, 2, 2);
    }

    return out;
}

/// Computes the conformal christoffel symbol - using tensors
AMREX_GPU_DEVICE AMREX_FORCE_INLINE chris_t compute_christoffel(
    const Tensor<3, amrex::Real> &d1_h, const Tensor<2, amrex::Real> &h_UU)
{
    chris_t out{};

    FOR (i, j, k)
    {
        out.LLL[i][j][k] =
            0.5 * (d1_h[j][i][k] + d1_h[k][i][j] - d1_h[j][k][i]);
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

/// Computes the conformal christoffel symbol - using tensors
AMREX_GPU_DEVICE AMREX_FORCE_INLINE chris_array_t
compute_christoffel(const amrex::Array3D<amrex::Real, 0, 3, 0, 3, 0, 3> &d1_h,
                    const amrex::Array2D<amrex::Real, 0, 3, 0, 3> &h_UU)
{
    chris_t_array out{};

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

/// Computes the conformal christoffel symbol - using tensors
AMREX_GPU_DEVICE AMREX_FORCE_INLINE chris_t_array
compute_christoffel(const amrex::Array3D<amrex::Real, 0, 3, 0, 3, 0, 3> &d1_h,
                    const amrex::Array1D<amrex::Real, 0, 6> &h_UU)
{
    chris_t_array out{};

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
            out.ULL(i, j, k) += h_UU(SYMM_IDX(i, l)) * out.LLL(l, j, k);
        }
    }
    FOR (i)
    {
        out.contracted(i) = 0;
        FOR (j, k)
        {
            out.contracted(i) += h_UU(SYMM_IDX(j, k)) * out.ULL(i, j, k);
        }
    }

    return out;
}

/// Computes the conformal christoffel symbol - using tensors
AMREX_GPU_DEVICE AMREX_FORCE_INLINE chris_t_array
compute_christoffel_test(const Tensor<3, amrex::Real> &d1_h,
                         const amrex::Array2D<amrex::Real, 0, 3, 0, 3> &h_UU)
{
    chris_t_array out{};

    FOR (i, j, k)
    {
        out.LLL(i, j, k) =
            0.5 * (d1_h[j][i][k] + d1_h[k][i][j] - d1_h[j][k][i]);
    }
    FOR (i, j, k)
    {
        out.ULL(i, j, k) = 0;
        FOR (l)
        {
            int idx           = i + l + ((i * l != 0) ? 1 : 0);
            out.ULL(i, j, k) += h_UU(i, l) * out.LLL(l, j, k);
        }
    }
    FOR (i)
    {
        out.contracted(i) = 0;
        FOR (j, k)
        {
            int idx            = j + k + ((j * k != 0) ? 1 : 0);
            out.contracted(i) += h_UU(j, k) * out.ULL(i, j, k);
        }
    }

    return out;
}

/// Computes the conformal christoffel symbol - using tensors
AMREX_GPU_DEVICE AMREX_FORCE_INLINE chris_t_array compute_christoffel_test(
    const amrex::Array3D<amrex::Real, 0, 3, 0, 3, 0, 3> &d1_h,
    const amrex::Array2D<amrex::Real, 0, 3, 0, 3> &h_UU)
{
    chris_t_array out{};

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
            int idx           = i + l + ((i * l != 0) ? 1 : 0);
            out.ULL(i, j, k) += h_UU(i, l) * out.LLL(l, j, k);
        }
    }
    FOR (i)
    {
        out.contracted(i) = 0;
        FOR (j, k)
        {
            int idx            = j + k + ((j * k != 0) ? 1 : 0);
            out.contracted(i) += h_UU(j, k) * out.ULL(i, j, k);
        }
    }

    return out;
}

/// Computes the conformal christoffel symbol
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE amrex::Array3D<amrex::Real, 0, 3, 0, 3, 0, 3>
compute_phys_chris(
    const Tensor<1, amrex::Real> &d1_chi, const CCZ4Vars &vars,
    const amrex::Array2D<amrex::Real, 0, 3, 0, 3> &h_UU,
    const amrex::Array3D<amrex::Real, 0, 3, 0, 3, 0, 3> &chris_ULL)
{
    using namespace TensorAlgebra;
    amrex::Array3D<amrex::Real, 0, 3, 0, 3, 0, 3> chris_phys;
    FOR (i, j, k)
    {
        chris_phys(i, j, k) =
            chris_ULL(i, j, k) -
            0.5 / vars.chi() *
                (delta(i, k) * d1_chi[j] + delta(i, j) * d1_chi[k]);
        FOR (m)
        {
            chris_phys(i, j, k) +=
                0.5 / vars.chi() * vars.h(j, k) * h_UU(i, m) * d1_chi[m];
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

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_divshift(const Tensor<2, amrex::Real> &d1_shift)
{
    amrex::Real divshift = 0.;
    FOR (i)
        divshift += d1_shift[i][i];
    return divshift;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_divshift(amrex::Array2D<amrex::Real, 0, 3, 0, 3> &d1_shift)
{
    amrex::Real divshift = 0.;
    FOR (i)
        divshift += d1_shift(i, i);
    return divshift;
}

/// Removes the trace of a 2-Tensor with lower indices given a metric and an
/// inverse metric.

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
make_trace_free(amrex::Array2D<amrex::Real, 0, 3, 0, 3> &tensor_LL,
                const CCZ4Vars vars,
                const amrex::Array2D<amrex::Real, 0, 3, 0, 3> &inverse_metric)
{
    auto trace = TensorAlgebra::compute_trace(tensor_LL, inverse_metric);
    double one_over_gr_spacedim = 1. / ((double)GR_SPACEDIM);
    FOR (i, j)
    {
        tensor_LL(i, j) -= one_over_gr_spacedim * vars.h(i, j) * trace;
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
make_trace_free(Tensor<2, amrex::Real> &tensor_LL,
                const amrex::CellData<const amrex::Real> h,
                const Tensor<2, amrex::Real> &inverse_metric)
{
    auto trace = TensorAlgebra::compute_trace(tensor_LL, inverse_metric);
    double one_over_gr_spacedim = 1. / ((double)GR_SPACEDIM);
    FOR (i, j)
    {
        tensor_LL[i][j] -=
            one_over_gr_spacedim * h[var_idx(c_h11, i, j)] * trace;
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
make_trace_free(amrex::Array2D<amrex::Real, 0, 3, 0, 3> &tensor_LL,
                const amrex::CellData<const amrex::Real> h,
                const amrex::Array1D<amrex::Real, 0, 6> &inverse_metric_sym)
{
    auto trace = TensorAlgebra::compute_trace(tensor_LL, inverse_metric_sym);
    double one_over_gr_spacedim = 1. / ((double)GR_SPACEDIM);
    FOR (i, j)
    {
        tensor_LL(i, j) -=
            one_over_gr_spacedim * h[var_idx(c_h11, i, j)] * trace;
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
make_trace_free(amrex::Array2D<amrex::Real, 0, 3, 0, 3> &tensor_LL,
                const amrex::CellData<const amrex::Real> h,
                const amrex::Array2D<amrex::Real, 0, 3, 0, 3> &inverse_metric)
{
    auto trace = TensorAlgebra::compute_trace(tensor_LL, inverse_metric);
    double one_over_gr_spacedim = 1. / ((double)GR_SPACEDIM);
    FOR (i, j)
    {
        tensor_LL(i, j) -=
            one_over_gr_spacedim * h[var_idx(c_h11, i, j)] * trace;
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
make_trace_free(Tensor<2, amrex::Real> &tensor_LL,
                const amrex::CellData<const amrex::Real> h,
                const amrex::Array2D<amrex::Real, 0, 3, 0, 3> &inverse_metric)
{
    auto trace = TensorAlgebra::compute_trace(tensor_LL, inverse_metric);
    double one_over_gr_spacedim = 1. / ((double)GR_SPACEDIM);
    FOR (i, j)
    {
        tensor_LL[i][j] -=
            one_over_gr_spacedim * h[var_idx(c_h11, i, j)] * trace;
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_z_terms(const int i, const int j,
                const amrex::Array1D<amrex::Real, 0, 3> &Z_over_chi,
                const CCZ4Vars &vars, const Tensor<1, amrex::Real> &d1_chi)
{
    amrex::Real out = 0.;
    FOR (k)
    {
        out += Z_over_chi(k) *
               (vars.h(i, k) * d1_chi[j] + vars.h(j, k) * d1_chi[i] -
                vars.h(i, j) * d1_chi[k]);
    }
    return out;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_z_terms(const int i, const int j,
                const Tensor<1, amrex::Real> &Z_over_chi,
                const amrex::CellData<const amrex::Real> &h,
                const Tensor<1, amrex::Real> &d1_chi)
{
    amrex::Real out = 0.;
    FOR (k)
    {
        out += Z_over_chi[k] * (h[var_idx(c_h11, i, k)] * d1_chi[j] +
                                h[var_idx(c_h11, j, k)] * d1_chi[i] -
                                h[var_idx(c_h11, i, j)] * d1_chi[k]);
    }
    return out;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_z_terms(const int i, const int j,
                const amrex::Array1D<amrex::Real, 0, 3> &Z_over_chi,
                const amrex::CellData<const amrex::Real> &h,
                const amrex::Array1D<amrex::Real, 0, 3> &d1_chi)
{
    amrex::Real out = 0.;
    FOR (k)
    {
        out += Z_over_chi(k) * (h[var_idx(c_h11, i, k)] * d1_chi(j) +
                                h[var_idx(c_h11, j, k)] * d1_chi(i) -
                                h[var_idx(c_h11, i, j)] * d1_chi(k));
    }

    return out;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
compute_z_terms(const int i, const int j,
                const Tensor<1, amrex::Real> &Z_over_chi,
                const amrex::CellData<const amrex::Real> &h,
                const amrex::Array1D<amrex::Real, 0, 3> &d1_chi)
{
    amrex::Real out = 0.;
    FOR (k)
    {
        out += Z_over_chi[k] * (h[var_idx(c_h11, i, k)] * d1_chi(j) +
                                h[var_idx(c_h11, j, k)] * d1_chi(i) -
                                h[var_idx(c_h11, i, j)] * d1_chi(k));
    }

    return out;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE ricci_t_array compute_ricci_Z(
    int ix, int iy, int iz, const amrex::Array4<amrex::Real> &rhs_state,
    const amrex::Array4<amrex::Real const> &state,
    const amrex::Array2D<amrex::Real, 0, 3, 0, 3> &h_UU,
    const chris_t_array &chris, const Tensor<1, amrex::Real> &Z_over_chi,
    const FourthOrderDerivatives &m_deriv)
{
    ricci_t_array out;

    const amrex::CellData<amrex::Real> &rhs_cell_data =
        rhs_state.cellData(ix, iy, iz);
    const amrex::CellData<const amrex::Real> &state_cell_data =
        state.cellData(ix, iy, iz);

    amrex::Array1D<amrex::Real, 0, 3> shift_vector;
    FOR (idir)
    {
        shift_vector(idir) = state(ix, iy, iz, c_shift1 + idir);
    }

    // chi derivatives
    Tensor<1, amrex::Real> d1_chi =
        m_deriv.diff1_scalar(ix, iy, iz, state, c_chi);

    auto d2_chi = m_deriv.diff2_sym_scalar(ix, iy, iz, state, c_chi);

    Tensor<2, amrex::Real> covdtilde2chi;
    FOR (k, l)
    {
        int idx             = k + l + ((k * l != 0) ? 1 : 0);
        covdtilde2chi[k][l] = d2_chi(idx);
        FOR (m)
        {
            covdtilde2chi[k][l] -= chris.ULL(m, k, l) * d1_chi[m];
        }
    }

    //    Tensor<3, amrex::Real> chris_LLU = {0.};
    amrex::Array3D<amrex::Real, 0, 3, 0, 3, 0, 3> chris_LLU;
    amrex::Real boxtildechi   = 0.;
    amrex::Real dchi_dot_dchi = 0.;
    FOR (i, j)
    {
        int idx        = i + j + ((i * j != 0) ? 1 : 0);
        boxtildechi   += covdtilde2chi[i][j] * h_UU(i, j);
        dchi_dot_dchi += d1_chi[i] * d1_chi[j] * h_UU(i, j);
    }

    // Gamma derivatives
    auto d1_Gamma = m_deriv.diff1_array_vector(ix, iy, iz, state, c_Gamma1);

    // hij derivatives
    auto d1_h = m_deriv.diff1_sym_tensor(ix, iy, iz, state, c_h11);

    // auto d2_h = m_deriv.diff2_sym_tensor(ix, iy, iz, state, c_h11);
    auto d2_h = m_deriv.diff2_sym_tensor_test_array(ix, iy, iz, state, c_h11);

    FOR (i, j)
    {
        amrex::Real ricci_hat = 0;
        int idx1              = i + j + ((i * j != 0) ? 1 : 0);

        FOR (k)
        {

            amrex::Real chris_LLU_jkl;
            amrex::Real chris_LLU_ikl;
            amrex::Real chris_LLU_kjl;

            // We call this ricci_hat rather than ricci_tilde as we have
            // replaced what should be \tilde{Gamma} with \hat{Gamma} in
            // order to avoid adding terms that cancel later on
            ricci_hat +=
                0.5 * (state_cell_data[var_idx(c_h11, k, i)] * d1_Gamma(k, j) +
                       state_cell_data[var_idx(c_h11, k, j)] * d1_Gamma(k, i));
            ricci_hat += 0.5 * state_cell_data[c_Gamma1 + k] * d1_h(idx1, k);

            FOR (l)
            {

                // jkl
                chris_LLU_jkl = h_UU(l, 0) * chris.LLL(j, k, 0) +
                                h_UU(l, 1) * chris.LLL(j, k, 1) +
                                h_UU(l, 2) * chris.LLL(j, k, 2);
                // ikl
                chris_LLU_ikl = h_UU(l, 0) * chris.LLL(i, k, 0) +
                                h_UU(l, 1) * chris.LLL(i, k, 1) +
                                h_UU(l, 2) * chris.LLL(i, k, 2);
                // kjl
                chris_LLU_kjl = h_UU(l, 0) * chris.LLL(k, j, 0) +
                                h_UU(l, 1) * chris.LLL(k, j, 1) +
                                h_UU(l, 2) * chris.LLL(k, j, 2);

                int idx2 = k + l + ((k * l != 0) ? 1 : 0);

                ricci_hat += -0.5 * h_UU(k, l) * d2_h(idx1, idx2) +
                             chris.ULL(k, l, i) * chris_LLU_jkl +
                             chris.ULL(k, l, j) * chris_LLU_ikl +
                             chris.ULL(k, i, l) * chris_LLU_kjl;
            }
        }

        amrex::Real ricci_chi =
            0.5 * ((GR_SPACEDIM - 2) * covdtilde2chi[i][j] +
                   state_cell_data[var_idx(c_h11, i, j)] * boxtildechi -
                   ((GR_SPACEDIM - 2) * d1_chi[i] * d1_chi[j] +
                    GR_SPACEDIM * state_cell_data[var_idx(c_h11, i, j)] *
                        dchi_dot_dchi) /
                       (2 * state_cell_data[c_chi]));

        amrex::Real z_terms =
            compute_z_terms(i, j, Z_over_chi, state_cell_data, d1_chi);

        out.LL(i, j) =
            (ricci_chi + state_cell_data[c_chi] * ricci_hat + z_terms) /
            state_cell_data[c_chi];
    }

    out.scalar =
        state_cell_data[c_chi] * TensorAlgebra::compute_trace(out.LL, h_UU);

    return out;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE ricci_t compute_ricci_Z_sym(
    int ix, int iy, int iz, const amrex::Array4<amrex::Real> &rhs_state,
    const amrex::Array4<amrex::Real const> &state,
    const amrex::Array1D<amrex::Real, 0, 6> &h_UU, const chris_t_array &chris,
    const Tensor<1, amrex::Real> &Z_over_chi,
    const FourthOrderDerivatives &m_deriv)
{
    ricci_t out;

    const amrex::CellData<amrex::Real> &rhs_cell_data =
        rhs_state.cellData(ix, iy, iz);
    const amrex::CellData<const amrex::Real> &state_cell_data =
        state.cellData(ix, iy, iz);

    // chi derivatives
    auto d1_chi = m_deriv.diff1_array_scalar(ix, iy, iz, state, c_chi);

    auto d2_chi = m_deriv.diff2_array_scalar(ix, iy, iz, state, c_chi);

    amrex::Array2D<amrex::Real, 0, 3, 0, 3> covdtilde2chi;
    FOR (k, l)
    {
        covdtilde2chi(k, l) = d2_chi(k, l);
        FOR (m)
        {
            covdtilde2chi(k, l) -= chris.ULL(m, k, l) * d1_chi(m);
        }
    }

    amrex::Array3D<amrex::Real, 0, 3, 0, 3, 0, 3> chris_LLU{};
    amrex::Real boxtildechi   = 0.;
    amrex::Real dchi_dot_dchi = 0.;
    FOR (i, j)
    {
        boxtildechi   += covdtilde2chi(i, j) * h_UU(SYMM_IDX(i, j));
        dchi_dot_dchi += d1_chi(i) * d1_chi(j) * h_UU(SYMM_IDX(i, j));
        FOR (k, l)
        {
            chris_LLU(i, j, k) += h_UU(SYMM_IDX(k, l)) * chris.LLL(i, j, l);
        }
    }

    // Gamma derivatives
    auto d1_Gamma = m_deriv.diff1_array_vector(ix, iy, iz, state, c_Gamma1);

    // hij derivatives
    auto d1_h = m_deriv.diff1_array_tensor(ix, iy, iz, state, c_h11);

    auto d2_h = m_deriv.diff2_sym_tensor(ix, iy, iz, state, c_h11);

    FOR (i, j)
    {
        amrex::Real ricci_hat = 0;
        int idx1              = i + j + ((i * j != 0) ? 1 : 0);
        FOR (k)
        {
            // We call this ricci_hat rather than ricci_tilde as we have
            // replaced what should be \tilde{Gamma} with \hat{Gamma} in
            // order to avoid adding terms that cancel later on

            ricci_hat +=
                0.5 * (state_cell_data[var_idx(c_h11, k, i)] * d1_Gamma(k, j) +
                       state_cell_data[var_idx(c_h11, k, j)] * d1_Gamma(k, i));
            ricci_hat += 0.5 * state_cell_data[c_Gamma1 + k] * d1_h(i, j, k);
            FOR (l)
            {
                int idx2   = k + l + ((k * l != 0) ? 1 : 0);
                ricci_hat += -0.5 * h_UU(SYMM_IDX(k, l)) * d2_h(i, j, idx2) +
                             (chris.ULL(k, l, i) * chris_LLU(j, k, l) +
                              chris.ULL(k, l, j) * chris_LLU(i, k, l) +
                              chris.ULL(k, i, l) * chris_LLU(k, j, l));
            }
        }

        amrex::Real ricci_chi =
            0.5 *
            ((GR_SPACEDIM - 2) * covdtilde2chi(i, j) +
             state_cell_data[c_h11 + idx1] * boxtildechi -
             ((GR_SPACEDIM - 2) * d1_chi(i) * d1_chi(j) +
              GR_SPACEDIM * state_cell_data[c_h11 + idx1] * dchi_dot_dchi) /
                 (2 * state_cell_data[c_chi]));

        amrex::Real z_terms =
            compute_z_terms(i, j, Z_over_chi, state_cell_data, d1_chi);

        out.LL[i][j] =
            (ricci_chi + state_cell_data[c_chi] * ricci_hat + z_terms) /
            state_cell_data[c_chi];
    }

    out.scalar =
        state_cell_data[c_chi] * TensorAlgebra::compute_trace(out.LL, h_UU);

    return out;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE ricci_t_array compute_ricci_Z(
    int ix, int iy, int iz, const amrex::Array4<amrex::Real> &rhs_state,
    const amrex::Array4<amrex::Real const> &state,
    const amrex::Array1D<amrex::Real, 0, 6> &h_UU, const chris_t_array &chris,
    const amrex::Array1D<amrex::Real, 0, 3> &Z_over_chi,
    const FourthOrderDerivatives &m_deriv)
{
    ricci_t_array out;

    const amrex::CellData<amrex::Real> &rhs_cell_data =
        rhs_state.cellData(ix, iy, iz);
    const amrex::CellData<const amrex::Real> &state_cell_data =
        state.cellData(ix, iy, iz);

    // amrex::Array1D<amrex::Real, 0, 3> shift_vector;
    // FOR (idir)
    // {
    //     shift_vector(idir) = state(ix, iy, iz, c_shift1 + idir);
    // }

    // chi derivatives
    auto d1_chi = m_deriv.diff1_array_scalar(ix, iy, iz, state, c_chi);

    auto d2_chi = m_deriv.diff2_sym_scalar(ix, iy, iz, state, c_chi);

    amrex::Array2D<amrex::Real, 0, 3, 0, 3> covdtilde2chi;
    FOR (k, l)
    {
        covdtilde2chi(k, l) = d2_chi(SYMM_IDX(k, l));
        FOR (m)
        {
            covdtilde2chi(k, l) -= chris.ULL(m, k, l) * d1_chi(m);
        }
    }

    amrex::Array3D<amrex::Real, 0, 3, 0, 3, 0, 3> chris_LLU{};
    amrex::Real boxtildechi   = 0.;
    amrex::Real dchi_dot_dchi = 0.;
    FOR (i, j)
    {
        boxtildechi   += covdtilde2chi(i, j) * h_UU(SYMM_IDX(i, j));
        dchi_dot_dchi += d1_chi(i) * d1_chi(j) * h_UU(SYMM_IDX(i, j));
        FOR (k, l)
        {
            chris_LLU(i, j, k) += h_UU(SYMM_IDX(k, l)) * chris.LLL(i, j, l);
        }
    }

    // Gamma derivatives
    auto d1_Gamma = m_deriv.diff1_array_vector(ix, iy, iz, state, c_Gamma1);

    // hij derivatives
    auto d1_h = m_deriv.diff1_sym_tensor(ix, iy, iz, state, c_h11);

    auto d2_h = m_deriv.diff2_sym_tensor(ix, iy, iz, state, c_h11);
    //    auto d2_h = m_deriv.diff2_sym_tensor_test_array(ix, iy, iz, state,
    //    c_h11);

    FOR (i, j)
    {
        amrex::Real ricci_hat = 0;
        int idx1              = i + j + ((i * j != 0) ? 1 : 0);
        FOR (k)
        {
            // We call this ricci_hat rather than ricci_tilde as we have
            // replaced what should be \tilde{Gamma} with \hat{Gamma} in
            // order to avoid adding terms that cancel later on

            ricci_hat +=
                0.5 * (state_cell_data[var_idx(c_h11, k, i)] * d1_Gamma(k, j) +
                       state_cell_data[var_idx(c_h11, k, j)] * d1_Gamma(k, i));
            ricci_hat += 0.5 * state_cell_data[c_Gamma1 + k] * d1_h(idx1, k);
            FOR (l)
            {

                int idx2   = k + l + ((k * l != 0) ? 1 : 0);
                ricci_hat += -0.5 * h_UU(idx2) * d2_h(i, j, idx2) +
                             (chris.ULL(k, l, i) * chris_LLU(j, k, l) +
                              chris.ULL(k, l, j) * chris_LLU(i, k, l) +
                              chris.ULL(k, i, l) * chris_LLU(k, j, l));
            }
        }

        amrex::Real ricci_chi =
            0.5 *
            ((GR_SPACEDIM - 2) * covdtilde2chi(i, j) +
             state_cell_data[c_h11 + idx1] * boxtildechi -
             ((GR_SPACEDIM - 2) * d1_chi(i) * d1_chi(j) +
              GR_SPACEDIM * state_cell_data[c_h11 + idx1] * dchi_dot_dchi) /
                 (2 * state_cell_data[c_chi]));

        amrex::Real z_terms =
            compute_z_terms(i, j, Z_over_chi, state_cell_data, d1_chi);

        out.LL(i, j) =
            (ricci_chi + state_cell_data[c_chi] * ricci_hat + z_terms) /
            state_cell_data[c_chi];
    }

    out.scalar =
        state_cell_data[c_chi] * TensorAlgebra::compute_trace(out.LL, h_UU);

    return out;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE ricci_t compute_ricci_Z(
    int ix, int iy, int iz, const amrex::Array4<amrex::Real> &rhs_state,
    const amrex::Array4<amrex::Real const> &state,
    const Tensor<2, amrex::Real> &h_UU, const chris_t &chris,
    const Tensor<1, amrex::Real> &Z_over_chi,
    const FourthOrderDerivatives &m_deriv)
{
    ricci_t out;

    const amrex::CellData<amrex::Real> &rhs_cell_data =
        rhs_state.cellData(ix, iy, iz);
    const amrex::CellData<const amrex::Real> &state_cell_data =
        state.cellData(ix, iy, iz);

    Tensor<1, amrex::Real> shift_vector;
    FOR (idir)
    {
        shift_vector[idir] = state(ix, iy, iz, c_shift1 + idir);
    }

    // chi derivatives
    Tensor<1, amrex::Real> d1_chi =
        m_deriv.diff1_scalar(ix, iy, iz, state, c_chi);

    auto d2_chi = m_deriv.diff2_scalar(ix, iy, iz, state, c_chi);

    Tensor<2, amrex::Real> covdtilde2chi;
    FOR (k, l)
    {
        covdtilde2chi[k][l] = d2_chi[k][l];
        FOR (m)
        {
            covdtilde2chi[k][l] -= chris.ULL[m][k][l] * d1_chi[m];
        }
    }

    Tensor<3, amrex::Real> chris_LLU = {0.};
    amrex::Real boxtildechi          = 0.;
    amrex::Real dchi_dot_dchi        = 0.;
    FOR (i, j)
    {
        int idx        = i + j + ((i * j != 0) ? 1 : 0);
        boxtildechi   += covdtilde2chi[i][j] * h_UU[i][j];
        dchi_dot_dchi += d1_chi[i] * d1_chi[j] * h_UU[i][j];
        FOR (k, l)
        {
            int idx             = k + l + ((k * l != 0) ? 1 : 0);
            chris_LLU[i][j][k] += h_UU[k][l] * chris.LLL[i][j][l];
        }
    }

    // Gamma derivatives
    Tensor<2, amrex::Real> d1_Gamma =
        m_deriv.diff1_vector(ix, iy, iz, state, c_Gamma1);

    // hij derivatives
    Tensor<3, amrex::Real> d1_h =
        m_deriv.diff1_tensor(ix, iy, iz, state, c_h11);

    // Tensor<4, amrex::Real> d2_h =
    //        m_deriv.diff2_tensor(ix, iy, iz, state, c_h11);
    Tensor<2, amrex::Real, 6> d2_h =
        m_deriv.diff2_sym_tensor_test(ix, iy, iz, state, c_h11);

    FOR (i, j)
    {
        amrex::Real ricci_hat = 0;
        int idx1              = i + j + ((i * j != 0) ? 1 : 0);
        FOR (k)
        {
            // We call this ricci_hat rather than ricci_tilde as we have
            // replaced what should be \tilde{Gamma} with \hat{Gamma} in
            // order to avoid adding terms that cancel later on

            ricci_hat +=
                0.5 * (state_cell_data[var_idx(c_h11, k, i)] * d1_Gamma[k][j] +
                       state_cell_data[var_idx(c_h11, k, j)] * d1_Gamma[k][i]);
            ricci_hat += 0.5 * state_cell_data[c_Gamma1 + k] * d1_h[i][j][k];
            FOR (l)
            {
                int idx2   = k + l + ((k * l != 0) ? 1 : 0);
                ricci_hat += -0.5 * h_UU[k][l] * d2_h[idx1][idx2] +
                             (chris.ULL[k][l][i] * chris_LLU[j][k][l] +
                              chris.ULL[k][l][j] * chris_LLU[i][k][l] +
                              chris.ULL[k][i][l] * chris_LLU[k][j][l]);
            }
        }

        amrex::Real ricci_chi =
            0.5 * ((GR_SPACEDIM - 2) * covdtilde2chi[i][j] +
                   state_cell_data[var_idx(c_h11, i, j)] * boxtildechi -
                   ((GR_SPACEDIM - 2) * d1_chi[i] * d1_chi[j] +
                    GR_SPACEDIM * state_cell_data[var_idx(c_h11, i, j)] *
                        dchi_dot_dchi) /
                       (2 * state_cell_data[c_chi]));

        amrex::Real z_terms =
            compute_z_terms(i, j, Z_over_chi, state_cell_data, d1_chi);

        out.LL[i][j] =
            (ricci_chi + state_cell_data[c_chi] * ricci_hat + z_terms) /
            state_cell_data[c_chi];
    }

    out.scalar =
        state_cell_data[c_chi] * TensorAlgebra::compute_trace(out.LL, h_UU);

    return out;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE ricci_t compute_ricci_Z(
    const CCZ4Vars &vars, const CCZ4D1Vars &d1,
    const Tensor<2, amrex::Real> &d2_chi, const Tensor<4, amrex::Real> &d2_h,
    const amrex::Array2D<amrex::Real, 0, 3, 0, 3> &h_UU,
    const chris_array_t &chris,
    const amrex::Array1D<amrex::Real, 0, 3> &Z_over_chi)
{
    ricci_t out;

    amrex::Array2D<amrex::Real, 0, 3, 0, 3> covdtilde2chi{};
    FOR (k, l)
    {
        covdtilde2chi(k, l) = d2_chi[k][l];
        FOR (m)
        {
            covdtilde2chi(k, l) -= chris.ULL(m, k, l) * d1.chi(m);
        }
    }

    amrex::Array3D<amrex::Real, 0, 3, 0, 3, 0, 3> chris_LLU{};
    amrex::Real boxtildechi   = 0.;
    amrex::Real dchi_dot_dchi = 0.;
    FOR (i, j)
    {
        boxtildechi   += covdtilde2chi(i, j) * h_UU(i, j);
        dchi_dot_dchi += d1.chi(i) * d1.chi(j) * h_UU(i, j);
        FOR (k, l)
        {
            chris_LLU(i, j, k) += h_UU(k, l) * chris.LLL(i, j, l);
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

                ricci_hat += -0.5 * h_UU(k, l) * d2_h[i][j][k][l] +
                             (chris.ULL(k, l, i) * chris_LLU(j, k, l) +
                              chris.ULL(k, l, j) * chris_LLU(i, k, l) +
                              chris.ULL(k, i, l) * chris_LLU(k, j, l));
            }
        }

        amrex::Real ricci_chi =
            0.5 * ((GR_SPACEDIM - 2) * covdtilde2chi(i, j) +
                   vars.h(i, j) * boxtildechi -
                   ((GR_SPACEDIM - 2) * d1.chi(i) * d1.chi(j) +
                    GR_SPACEDIM * vars.h(i, j) * dchi_dot_dchi) /
                       (2 * vars.chi()));

        amrex::Real z_terms = compute_z_terms(i, j, Z_over_chi, vars, d1.chi());

        out.LL(i, j) =
            (ricci_chi + vars.chi() * ricci_hat + z_terms) / vars.chi();
    }

    out.scalar = vars.chi() * TensorAlgebra::compute_trace(out.LL, h_UU);

    return out;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Array2D<amrex::Real, 0, 3, 0, 3>
compute_d1_chris_contracted(const amrex::Array2D<amrex::Real, 0, 3, 0, 3> &h_UU,
                            const CCZ4D1Vars &d1,
                            const Tensor<4, amrex::Real> &d2_h)
{
    amrex::Array2D<amrex::Real, 0, 3, 0, 3> d1_chris_contracted{};

    FOR (i, j)
    {
        FOR (m, n, p)
        {
            amrex::Real d1_terms = 0.0;
            FOR (q, r)
            {
                d1_terms += -h_UU(q, r) * (d1.h(n, q, j) * d1.h(m, p, r) +
                                           d1.h(m, n, j) * d1.h(p, q, r));
            }
            d1_chris_contracted(i, j) +=
                h_UU(i, m) * h_UU(n, p) * (d2_h[m][n][j][p] + d1_terms);
        }
    }
    return d1_chris_contracted;
}

// This function allows adding arbitrary multiples of D_{(i}Z_{j)}
// to the Ricci scalar rather than the default of 2 in compute_ricci_Z
AMREX_GPU_DEVICE AMREX_FORCE_INLINE ricci_t compute_ricci_Z_general(
    const CCZ4Vars &vars, const CCZ4D1Vars &d1,
    const Tensor<2, amrex::Real> &d2_chi, const Tensor<4, amrex::Real> &d2_h,
    const amrex::Array2D<amrex::Real, 0, 3, 0, 3> &h_UU,
    const chris_array_t &chris, const double dZ_coeff)
{
    // get contributions from conformal metric and factor with zero Z vector
    amrex::Array1D<amrex::Real, 0, 3> zero_Z{};
    auto ricci = compute_ricci_Z(vars, d1, d2_chi, d2_h, h_UU, chris, zero_Z);

    // need to add term to correct for d1.Gamma (includes Z contribution)
    // and Gamma in ricci_hat
    auto d1_chris_contracted = compute_d1_chris_contracted(h_UU, d1, d2_h);
    amrex::Array1D<amrex::Real, 0, 3> Z_over_chi;
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
                (vars.h(m, i) * (d1_chris_contracted(m, j) - d1.Gamma(m, j)) +
                 vars.h(m, j) * (d1_chris_contracted(m, i) - d1.Gamma(m, i)) +
                 (chris.contracted(m) - vars.Gamma(m)) * d1.h(i, j, m));
        }
        amrex::Real z_terms = compute_z_terms(i, j, Z_over_chi, vars, d1.chi());
        ricci.LL(i, j) += 0.5 * dZ_coeff * z_terms / vars.chi();
    }
    ricci.scalar = vars.chi() * TensorAlgebra::compute_trace(ricci.LL, h_UU);
    return ricci;
}

// This function returns the pure Ricci scalar with no contribution from the
// Z vector - used e.g. in the constraint calculations.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE ricci_t compute_ricci(
    const CCZ4Vars &vars, const CCZ4D1Vars &d1,
    const Tensor<2, amrex::Real> &d2_chi, const Tensor<4, amrex::Real> &d2_h,
    const amrex::Array2D<amrex::Real, 0, 3, 0, 3> &h_UU,
    const chris_array_t &chris)
{
    return compute_ricci_Z_general(vars, d1, d2_chi, d2_h, h_UU, chris, 0.0);
}

} // namespace CCZ4Geometry

#endif /* CCZ4GEOMETRY_HPP_ */
