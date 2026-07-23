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

[[nodiscard]] constexpr inline int sym_var_idx(const int ivar, const int i,
                                               const int j) noexcept
{
    return ivar + i + j + ((i * j != 0) ? 1 : 0);
}

[[nodiscard]] constexpr inline int sym_var_idx(const int i,
                                               const int j) noexcept
{
    return i + j + ((i * j != 0) ? 1 : 0);
}

#define SPACETIME_DIM GR_SPACEDIM + 1

// Number of unique indices after accounting for symmetry
#define NUM_SYM_IDXS (AMREX_SPACEDIM * (AMREX_SPACEDIM + 1) / 2)

namespace Tensor
{

template <int rank, int... DIMS>
requires(sizeof...(DIMS) == rank && rank <= 4) AMREX_GPU_HOST_DEVICE
    struct GeneralRank;

template <int DIM> AMREX_GPU_HOST_DEVICE struct GeneralRank<1, DIM>
{
  private:
    amrex::Array1D<amrex::Real, 0, DIM - 1> m_tensor{};

  public:
    AMREX_GPU_HOST_DEVICE GeneralRank() = default;

    AMREX_GPU_HOST_DEVICE
    GeneralRank(std::initializer_list<amrex::Real> some_initial_values)
    {
        for (int i = 0; i < some_initial_values.size(); ++i)
        {
            m_tensor(i) = some_initial_values.begin()[i];
        }
    }

    AMREX_GPU_HOST_DEVICE amrex::Real &operator()(int idx)
    {
        return m_tensor(idx);
    }
    AMREX_GPU_HOST_DEVICE const amrex::Real &operator()(int idx) const
    {
        return m_tensor(idx);
    }

    AMREX_GPU_HOST_DEVICE amrex::Real &operator()(int idx1, int idx2)
        requires(DIM == NUM_SYM_IDXS)
    {
        return m_tensor(sym_var_idx(idx1, idx2));
    }

    AMREX_GPU_HOST_DEVICE const amrex::Real &operator()(int idx1,
                                                        int idx2) const
        requires(DIM == NUM_SYM_IDXS)
    {
        return m_tensor(sym_var_idx(idx1, idx2));
    }

    AMREX_GPU_HOST_DEVICE GeneralRank<1, DIM> &
    operator=(const amrex::Real an_initial_value)
    {
        for (int i = 0; i < DIM; ++i)
        {
            m_tensor(i) = an_initial_value;
        }

        return *this;
    }

    AMREX_GPU_HOST_DEVICE GeneralRank<1, DIM> &
    operator=(std::initializer_list<amrex::Real> some_initial_values)
    {
        for (int i = 0; i < some_initial_values.size(); ++i)
        {
            m_tensor(i) = some_initial_values.begin()[i];
        }

        return *this;
    }
};

template <int DIM1, int DIM2>
AMREX_GPU_HOST_DEVICE struct GeneralRank<2, DIM1, DIM2>
{
  private:
    amrex::Array2D<amrex::Real, 0, DIM1 - 1, 0, DIM2 - 1> m_tensor{};

  public:
    AMREX_GPU_HOST_DEVICE GeneralRank() = default;

    AMREX_GPU_HOST_DEVICE GeneralRank(amrex::Real an_initial_value)
    {
        for (int i = 0; i < DIM1; ++i)
            for (int j = 0; j < DIM2; ++j)
            {
                m_tensor(i, j) = an_initial_value;
            }
    }

    AMREX_GPU_HOST_DEVICE amrex::Real &operator()(int idx1, int idx2)
    {
        return m_tensor(idx1, idx2);
    }

    AMREX_GPU_HOST_DEVICE const amrex::Real &operator()(int idx1,
                                                        int idx2) const
    {
        return m_tensor(idx1, idx2);
    }

    AMREX_GPU_HOST_DEVICE amrex::Real &operator()(int idx1, int idx2, int idx3,
                                                  int idx4)
        requires(DIM1 == NUM_SYM_IDXS && DIM2 == NUM_SYM_IDXS)
    {
        return m_tensor(sym_var_idx(idx1, idx2), sym_var_idx(idx3, idx4));
    }

    AMREX_GPU_HOST_DEVICE const amrex::Real &
    operator()(int idx1, int idx2, int idx3, int idx4) const
        requires(DIM1 == NUM_SYM_IDXS && DIM2 == NUM_SYM_IDXS)
    {
        return m_tensor(sym_var_idx(idx1, idx2), sym_var_idx(idx3, idx4));
    }

    AMREX_GPU_HOST_DEVICE amrex::Real &operator()(int idx1, int idx2, int idx3)
        requires(DIM1 == NUM_SYM_IDXS || DIM2 == NUM_SYM_IDXS)
    {

        if constexpr (DIM1 == NUM_SYM_IDXS)
        {
            return m_tensor(sym_var_idx(idx1, idx2), idx3);
        }
        else
        {
            return m_tensor(idx1, sym_var_idx(idx2, idx3));
        }
    }

    AMREX_GPU_HOST_DEVICE const amrex::Real &operator()(int idx1, int idx2,
                                                        int idx3) const
        requires(DIM1 == NUM_SYM_IDXS || DIM2 == NUM_SYM_IDXS)
    {
        if constexpr (DIM1 == NUM_SYM_IDXS)
        {
            return m_tensor(sym_var_idx(idx1, idx2), idx3);
        }
        else
        {
            return m_tensor(idx1, sym_var_idx(idx2, idx3));
        }
    }

    AMREX_GPU_HOST_DEVICE GeneralRank<2, DIM1, DIM2> &
    operator=(const amrex::Real a_value)
    {
        for (int i = 0; i < DIM1; ++i)
            for (int j = 0; j < DIM2; ++j)
            {
                m_tensor(i, j) = a_value;
            }
        return *this;
    }
};

template <int DIM1, int DIM2, int DIM3>
AMREX_GPU_HOST_DEVICE struct GeneralRank<3, DIM1, DIM2, DIM3>
{
  private:
    amrex::Array3D<amrex::Real, 0, DIM1 - 1, 0, DIM2 - 1, 0, DIM3 - 1>
        m_tensor{};

  public:
    AMREX_GPU_HOST_DEVICE GeneralRank() = default;

    AMREX_GPU_HOST_DEVICE GeneralRank(amrex::Real an_initial_value)
    {
        for (int i = 0; i < DIM1; ++i)
            for (int j = 0; j < DIM2; ++j)
                for (int k = 0; k < DIM3; ++k)
                {
                    m_tensor(i, j, k) = an_initial_value;
                }
    }

    AMREX_GPU_HOST_DEVICE amrex::Real &operator()(int idx1, int idx2, int idx3)
    {
        return m_tensor(idx1, idx2, idx3);
    }
    AMREX_GPU_HOST_DEVICE const amrex::Real &operator()(int idx1, int idx2,
                                                        int idx3) const
    {
        return m_tensor(idx1, idx2, idx3);
    }

    AMREX_GPU_HOST_DEVICE GeneralRank<3, DIM1, DIM2, DIM3> &
    operator=(const amrex::Real a_value)
    {
        for (int i = 0; i < DIM1; ++i)
            for (int j = 0; j < DIM2; ++j)
                for (int k = 0; j < DIM3; ++k)
                {
                    m_tensor(i, j, k) = a_value;
                }

        return *this;
    }
};

template <int DIM1, int DIM2, int DIM3, int DIM4>
AMREX_GPU_HOST_DEVICE struct GeneralRank<4, DIM1, DIM2, DIM3, DIM4>
{
  private:
    amrex::Array1D<amrex::Real, 0, DIM1 * DIM2 * DIM3 * DIM4 - 1> m_tensor{};

  public:
    AMREX_GPU_HOST_DEVICE GeneralRank() = default;

    AMREX_GPU_HOST_DEVICE GeneralRank(amrex::Real an_initial_value)
    {
        const int total_dims = DIM1 * DIM2 * DIM3 * DIM4 - 1;

        for (int i = 0; i < total_dims; ++i)
            m_tensor(i) = an_initial_value;
    }

    AMREX_GPU_HOST_DEVICE amrex::Real &operator()(const int idx)
    {
        return m_tensor(idx);
    }

    AMREX_GPU_HOST_DEVICE const amrex::Real &operator()(const int idx) const
    {
        return m_tensor(idx);
    }

    AMREX_GPU_HOST_DEVICE amrex::Real &operator()(int idx1, int idx2, int idx3,
                                                  int idx4)
    {
        const int idx =
            DIM1 * DIM2 * DIM3 * idx1 + DIM2 * DIM3 * idx2 + DIM3 * idx3 + idx4;
        return m_tensor(idx);
    }
    AMREX_GPU_HOST_DEVICE const amrex::Real &
    operator()(int idx1, int idx2, int idx3, int idx4) const
    {
        const int idx =
            DIM1 * DIM2 * DIM3 * idx1 + DIM2 * DIM3 * idx2 + DIM3 * idx3 + idx4;
        return m_tensor(idx);
    }

    AMREX_GPU_HOST_DEVICE
    GeneralRank<4, DIM1, DIM2, DIM3, DIM4> &
    operator=(const amrex::Real &a_value)
    {
        const int total_dims = DIM1 * DIM2 * DIM3 * DIM4 - 1;

        for (int i = 0; i < total_dims; ++i)
            m_tensor(i) = a_value;

        return *this;
    }
};

// The default dimension use AMREX_SPACEDIM aka TENSOR_DIM if not otherwise
// defined
using Rank1 = GeneralRank<1, AMREX_SPACEDIM>;
using Rank2 = GeneralRank<2, AMREX_SPACEDIM, AMREX_SPACEDIM>;
using Rank3 = GeneralRank<3, AMREX_SPACEDIM, AMREX_SPACEDIM, AMREX_SPACEDIM>;
using Rank4 = GeneralRank<4, AMREX_SPACEDIM, AMREX_SPACEDIM, AMREX_SPACEDIM,
                          AMREX_SPACEDIM>;

// These are for symmetric tensors
using Sym12Rank2      = GeneralRank<1, NUM_SYM_IDXS>;
using Sym12Sym34Rank4 = GeneralRank<2, NUM_SYM_IDXS, NUM_SYM_IDXS>;
using Sym12Rank3      = GeneralRank<2, NUM_SYM_IDXS, AMREX_SPACEDIM>;
using Sym23Rank3      = GeneralRank<2, AMREX_SPACEDIM, NUM_SYM_IDXS>;

using SpaceTime =
    GeneralRank<4, SPACETIME_DIM, SPACETIME_DIM, SPACETIME_DIM, SPACETIME_DIM>;
} // namespace Tensor

#endif /* TENSOR_HPP_ */
