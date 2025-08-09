/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef TENSORALGEBRA2_HPP_
#define TENSORALGEBRA2_HPP_

#include "CCZ4Vars2.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include <AMReX_REAL.H>
#include <array>

namespace TensorAlgebra2
{
/// Computes the (i,j) component of the Kronecker delta
constexpr int delta(int i, int j) { return static_cast<int>(i == j); }

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
compute_trace_A(const CCZ4Vars2 &vars,
                const Tensor<2, amrex::Real> &inverse_metric)
{
    amrex::Real trace = 0.0;
    FOR (i, j)
    {
        trace += inverse_metric[i][j] * vars.A(i, j);
    }
    return trace;
}

} // namespace TensorAlgebra2

#endif /* TENSORALGEBRA2_HPP_ */
