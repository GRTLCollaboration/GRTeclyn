/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef DERIVATIVEBASE_HPP_
#define DERIVATIVEBASE_HPP_

#include "DimensionDefinitions.hpp"

#include <AMReX_REAL.H>
#include "AMReX_Array.H"

class DerivativeBase
{
  protected:
    amrex::Real m_dx;
    amrex::Real m_one_over_dx;
    amrex::Real m_one_over_dx2;

    AMREX_GPU_HOST_DEVICE
    DerivativeBase(amrex::Real dx)
        : m_dx(dx), m_one_over_dx(1 / dx), m_one_over_dx2(1 / (dx * dx))
    {
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE static const amrex::Real *
    get_var_ptr(const int ivar, const amrex::Real *state_ptr_xyz,
                const amrex::GpuArray<int, AMREX_SPACEDIM + 1> strides) noexcept
    {
        return state_ptr_xyz + ivar * strides[3];
    }

    [[nodiscard]] AMREX_GPU_DEVICE
        AMREX_FORCE_INLINE static amrex::GpuArray<int, AMREX_SPACEDIM + 1>
        get_strides(const amrex::Array4<const amrex::Real> &state) noexcept
    {
        int j_stride = static_cast<int>(state.stride.a[0]);
        int k_stride = static_cast<int>(state.stride.a[1]);
        int n_stride = static_cast<int>(state.stride.a[2]);

        amrex::GpuArray<int, AMREX_SPACEDIM + 1> strides{1, j_stride, k_stride,
                                                         n_stride};

        return strides;
    }
};

#endif /* DERIVATIVEBASE_HPP_ */