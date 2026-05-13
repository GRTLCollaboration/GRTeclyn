/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef TENSOR_HPP_
#define TENSOR_HPP_

#include "AlwaysInline.hpp"
#include "DimensionDefinitions.hpp"

// AMReX includes
#include "AMReX_Array.H"
#include "AMReX_GpuQualifiers.H"

// A function to return the right index for the tensors based on the
// ordering below 0: T11, 1: T12, 2: T13, 3: T22, 4: T23, 5: T33
#define VAR_IDX(ivar, i, j) ((ivar) + (i) + (j) + (((i) * (j) != 0) ? 1 : 0))

// A version for where the base reference for the tensor is 0
#define VAR_IDX0(i, j) VAR_IDX(0, (i), (j))

#define SPACETIME_DIM GR_SPACEDIM + 1

// Number of unique indices after accounting for symmetry
#define UNIQUE_IDX 6

namespace TensorArray
{
using Rank1 = amrex::Array1D<amrex::Real, 0, AMREX_SPACEDIM - 1>;
using Rank2 =
    amrex::Array2D<amrex::Real, 0, AMREX_SPACEDIM - 1, 0, AMREX_SPACEDIM - 1>;
using Rank3 = amrex::Array3D<amrex::Real, 0, AMREX_SPACEDIM - 1, 0,
                             AMREX_SPACEDIM - 1, 0, AMREX_SPACEDIM - 1>;

using Rank1Sym = amrex::Array1D<amrex::Real, 0, UNIQUE_IDX - 1>;
using Rank2Sym =
    amrex::Array2D<amrex::Real, 0, UNIQUE_IDX - 1, 0, UNIQUE_IDX - 1>;

// Rank 4 tensors are actually a 1D array of length 4 * 4 * 4 * 4

using Rank4 = amrex::Array1D<amrex::Real, 0, 255>;
} // namespace TensorArray

namespace SpaceTimeTensor
{
using Rank1 = amrex::Array1D<amrex::Real, 0, SPACETIME_DIM - 1>;
using Rank2 =
    amrex::Array2D<amrex::Real, 0, SPACETIME_DIM - 1, 0, SPACETIME_DIM - 1>;
using Rank3 = amrex::Array3D<amrex::Real, 0, SPACETIME_DIM - 1, 0,
                             SPACETIME_DIM - 1, 0, SPACETIME_DIM - 1>;

} // namespace SpaceTimeTensor

#endif /* TENSOR_HPP_ */
