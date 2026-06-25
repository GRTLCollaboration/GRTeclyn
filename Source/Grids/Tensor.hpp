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
#define NUM_SYM_IDXS (AMREX_SPACEDIM * (AMREX_SPACEDIM + 1) / 2)

namespace TensorArray
{
using Rank1 = amrex::Array1D<amrex::Real, 0, AMREX_SPACEDIM - 1>;
using Rank2 =
    amrex::Array2D<amrex::Real, 0, AMREX_SPACEDIM - 1, 0, AMREX_SPACEDIM - 1>;
using Rank3 = amrex::Array3D<amrex::Real, 0, AMREX_SPACEDIM - 1, 0,
                             AMREX_SPACEDIM - 1, 0, AMREX_SPACEDIM - 1>;

using Rank1Sym = amrex::Array1D<amrex::Real, 0, NUM_SYM_IDXS - 1>;
using Rank2Sym =
    amrex::Array2D<amrex::Real, 0, NUM_SYM_IDXS - 1, 0, NUM_SYM_IDXS - 1>;

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

namespace Tensor
{
template <int DIM = AMREX_SPACEDIM> struct GeneralRank1
{
    amrex::Array1D<amrex::Real, 0, DIM - 1> m_tensor{};

    amrex::Real &operator()(int idx) { return m_tensor(idx); }
    const amrex::Real &operator()(int idx) const { return m_tensor(idx); }
};

template <int DIM1 = AMREX_SPACEDIM, int DIM2 = AMREX_SPACEDIM>
struct GeneralRank2
{
    amrex::Array2D<amrex::Real, 0, DIM1 - 1, 0, DIM2 - 1> m_tensor{};

    amrex::Real &operator()(int idx1, int idx2) { return m_tensor(idx1, idx2); }
    const amrex::Real &operator()(int idx1, int idx2) const
    {
        return m_tensor(idx1, idx2);
    }
};

template <int DIM1 = AMREX_SPACEDIM, int DIM2 = AMREX_SPACEDIM,
          int DIM3 = AMREX_SPACEDIM>
struct GeneralRank3
{
    amrex::Array3D<amrex::Real, 0, DIM1 - 1, 0, DIM2 - 1, 0, DIM3 - 1>
        m_tensor{};

    amrex::Real &operator()(int idx1, int idx2, int idx3)
    {
        return m_tensor(idx1, idx2, idx3);
    }
    const amrex::Real &operator()(int idx1, int idx2, int idx3) const
    {
        return m_tensor(idx1, idx2, idx3);
    }
};

template <int DIM1 = SPACETIME_DIM, int DIM2 = SPACETIME_DIM,
          int DIM3 = SPACETIME_DIM, int DIM4 = SPACETIME_DIM>
struct GeneralRank4
{
    amrex::Array1D<amrex::Real, 0, DIM1 * DIM2 * DIM3 * DIM4 - 1> m_tensor{};

    amrex::Real &operator()(int idx1, int idx2, int idx3, int idx4)
    {
        const int idx = DIM1 * DIM2 * DIM3 * idx1 + DIM2 * DIM3 * idx2 +
                        DIM3 * idx3 + idx4; // TODO: check dimensions!
        return m_tensor(idx);
    }
    const amrex::Real &operator()(int idx1, int idx2, int idx3, int idx4) const
    {
        const int idx = DIM1 * DIM2 * DIM3 * idx1 + DIM2 * DIM3 * idx2 +
                        DIM3 * idx3 + idx4; // TODO: check dimensions!
        return m_tensor(idx);
    }
};

// The default dimension use AMREX_SPACEDIM aka TENSOR_DIM if not otherwise
// defined
using Rank1 = GeneralRank1<AMREX_SPACEDIM>;
using Rank2 = GeneralRank2<AMREX_SPACEDIM, AMREX_SPACEDIM>;
using Rank3 = GeneralRank3<AMREX_SPACEDIM, AMREX_SPACEDIM, AMREX_SPACEDIM>;
using Rank4 = GeneralRank4<AMREX_SPACEDIM, AMREX_SPACEDIM, AMREX_SPACEDIM,
                           AMREX_SPACEDIM>;

// These are for symmetric tensors
using Sym12Rank2      = GeneralRank1<NUM_SYM_IDXS>;
using Sym12Sym34Rank4 = GeneralRank2<NUM_SYM_IDXS, NUM_SYM_IDXS>;
using Sym12Rank3      = GeneralRank2<NUM_SYM_IDXS, AMREX_SPACEDIM>;
using Sym23Rank3      = GeneralRank2<AMREX_SPACEDIM, NUM_SYM_IDXS>;

using SpaceTime =
    GeneralRank4<SPACETIME_DIM, SPACETIME_DIM, SPACETIME_DIM, SPACETIME_DIM>;
} // namespace Tensor

template <typename T> class TensorAlt
{
  private:
    const int m_length{};
    T *m_tensor{};

  public:
    template <T... lengths>
    TensorAlt()
        : m_length((1 * ... * lengths)),
          m_tensor(amrex::Array1D<amrex::Real, 0, m_length>{}){};

    amrex::Real &operator()(int args, ...)
    {
        int index{1};
        //        index = lengths*index + lengths*args;
        // index = A0;  //1D
        // index = L1*A0 + A1; //2D
        // index = L1*L2*A1 + L2 * A1 + A3; //3D
        return m_tensor(index);
    }
};
#endif /* TENSOR_HPP_ */
