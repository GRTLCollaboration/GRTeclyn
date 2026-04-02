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

namespace TensorArray
{
using Rank1 = amrex::Array1D<amrex::Real, 0, 3>;
using Rank2 = amrex::Array2D<amrex::Real, 0, 3, 0, 3>;
using Rank3 = amrex::Array3D<amrex::Real, 0, 3, 0, 3, 0, 3>;

using Rank1Sym = amrex::Array1D<amrex::Real, 0, 6>;
using Rank2Sym = amrex::Array2D<amrex::Real, 0, 6, 0, 6>;

using Rank4 = amrex::Array1D<amrex::Real, 0, 256>;
} // namespace TensorArray

#endif /* TENSOR_HPP_ */
