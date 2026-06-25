/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef FOURTHORDERDERIVATIVES_HPP_
#define FOURTHORDERDERIVATIVES_HPP_

#include "Cell.hpp"
#include "DimensionDefinitions.hpp"

#include "StateVariables.hpp"
#include "Tensor.hpp"
#include <AMReX_REAL.H>
#include <array>
#include "AMReX_Array.H"

using namespace amrex::literals;

class FourthOrderDerivatives
{
  private:
    amrex::Real m_dx;
    amrex::Real m_one_over_dx;
    amrex::Real m_one_over_dx2;

  public:
    AMREX_GPU_HOST_DEVICE FourthOrderDerivatives(double dx)
        : m_dx(dx), m_one_over_dx(1 / dx), m_one_over_dx2(1 / (dx * dx))
    {
    }

    // NOLINTBEGIN(bugprone-easily-swappable-parameters)

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    diff1(const amrex::Real *in_ptr, const int stride, const int idx = 0) const
    {
        amrex::Real weight_far  = 8.33333333333333333333e-2_rt;
        amrex::Real weight_near = 6.66666666666666666667e-1_rt;

        return (weight_far * in_ptr[idx - 2 * stride] -
                weight_near * in_ptr[idx - stride] +
                weight_near * in_ptr[idx + stride] -
                weight_far * in_ptr[idx + 2 * stride]) *
               m_one_over_dx;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE TensorArray::Rank1
    diff1_scalar(int ix, int iy, int iz,
                 const amrex::Array4<const amrex::Real> &state,
                 const int ivar) const
    {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
        TensorArray::Rank1 d1;
        const auto *state_ptr_xyz = state.ptr(ix, iy, iz);
        amrex::GpuArray<int, AMREX_SPACEDIM> strides{
            1, static_cast<int>(state.stride.a[0]),
            static_cast<int>(state.stride.a[1])};
        const auto *var_ptr = state_ptr_xyz + ivar * state.stride.a[2];
        FOR (idir)
        {
            d1(idir) = diff1(var_ptr, strides[idir]);
        }
        return d1;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE TensorArray::Rank2
    diff1_vector(int ix, int iy, int iz,
                 const amrex::Array4<const amrex::Real> &state,
                 const int ivar_0) const
    {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
        TensorArray::Rank2 d1;
        const auto *state_ptr_xyz = state.ptr(ix, iy, iz);
        amrex::GpuArray<int, AMREX_SPACEDIM> strides{
            1, static_cast<int>(state.stride.a[0]),
            static_cast<int>(state.stride.a[1])};

        FOR (icomp)
        {
            const int ivar      = ivar_0 + icomp;
            const auto *var_ptr = state_ptr_xyz + ivar * state.stride.a[2];
            FOR (idir)
            {
                d1(icomp, idir) = diff1(var_ptr, strides[idir]);
            }
        }
        return d1;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        amrex::Array2D<amrex::Real, 0, NUM_SYM_IDXS - 1, 0, AMREX_SPACEDIM - 1>
        diff1_sym_tensor(int ix, int iy, int iz,
                         const amrex::Array4<const amrex::Real> &state,
                         const int ivar_0) const
    {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
        amrex::Array2D<amrex::Real, 0, NUM_SYM_IDXS - 1, 0, AMREX_SPACEDIM - 1>
            d1;
        const auto *state_ptr_xyz = state.ptr(ix, iy, iz);
        amrex::GpuArray<int, AMREX_SPACEDIM> strides{
            1, static_cast<int>(state.stride.a[0]),
            static_cast<int>(state.stride.a[1])};

        for (int ivar = 0; ivar < NUM_SYM_IDXS; ++ivar)
        {
            const auto *var_ptr =
                state_ptr_xyz + (ivar_0 + ivar) * state.stride.a[2];

            FOR (idir)
            {
                d1(ivar, idir) = diff1(var_ptr, strides[idir]);
            }
        }
        return d1;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE TensorArray::Rank3
    diff1_tensor(int ix, int iy, int iz,
                 const amrex::Array4<const amrex::Real> &state,
                 const int ivar_0) const
    {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
        TensorArray::Rank3 d1{};
        const auto *state_ptr_xyz = state.ptr(ix, iy, iz);
        amrex::GpuArray<int, AMREX_SPACEDIM> strides{
            1, static_cast<int>(state.stride.a[0]),
            static_cast<int>(state.stride.a[1])};

        FOR (icomp, jcomp)
        {
            const int ivar      = ivar_0 + VAR_IDX0(icomp, jcomp);
            const auto *var_ptr = state_ptr_xyz + ivar * state.stride.a[2];
            FOR (idir)
            {
                d1(icomp, jcomp, idir) = diff1(var_ptr, strides[idir]);
            }
        }
        return d1;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor::Rank1
    diff1_scalar_test(int ix, int iy, int iz,
                      const amrex::Array4<const amrex::Real> &state,
                      const int ivar) const
    {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
        Tensor::Rank1 d1;
        const auto *state_ptr_xyz = state.ptr(ix, iy, iz);
        amrex::GpuArray<int, AMREX_SPACEDIM> strides{
            1, static_cast<int>(state.stride.a[0]),
            static_cast<int>(state.stride.a[1])};
        const auto *var_ptr = state_ptr_xyz + ivar * state.stride.a[2];
        FOR (idir)
        {
            d1(idir) = diff1(var_ptr, strides[idir]);
        }
        return d1;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor::Rank2
    diff1_vector_test(int ix, int iy, int iz,
                      const amrex::Array4<const amrex::Real> &state,
                      const int ivar_0) const
    {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
        Tensor::Rank2 d1;
        const auto *state_ptr_xyz = state.ptr(ix, iy, iz);
        amrex::GpuArray<int, AMREX_SPACEDIM> strides{
            1, static_cast<int>(state.stride.a[0]),
            static_cast<int>(state.stride.a[1])};

        FOR (icomp)
        {
            const int ivar      = ivar_0 + icomp;
            const auto *var_ptr = state_ptr_xyz + ivar * state.stride.a[2];
            FOR (idir)
            {
                d1(icomp, idir) = diff1(var_ptr, strides[idir]);
            }
        }
        return d1;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor::Sym12Rank3
    diff1_sym_tensor_test(int ix, int iy, int iz,
                          const amrex::Array4<const amrex::Real> &state,
                          const int ivar_0) const
    {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
        Tensor::Sym12Rank3 d1;
        const auto *state_ptr_xyz = state.ptr(ix, iy, iz);
        amrex::GpuArray<int, AMREX_SPACEDIM> strides{
            1, static_cast<int>(state.stride.a[0]),
            static_cast<int>(state.stride.a[1])};

        for (int ivar = 0; ivar < NUM_SYM_IDXS; ++ivar)
        {
            const auto *var_ptr =
                state_ptr_xyz + (ivar_0 + ivar) * state.stride.a[2];

            FOR (idir)
            {
                d1(ivar, idir) = diff1(var_ptr, strides[idir]);
            }
        }
        return d1;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor::Rank3
    diff1_tensor_test(int ix, int iy, int iz,
                      const amrex::Array4<const amrex::Real> &state,
                      const int ivar_0) const
    {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
        Tensor::Rank3 d1{};
        const auto *state_ptr_xyz = state.ptr(ix, iy, iz);
        amrex::GpuArray<int, AMREX_SPACEDIM> strides{
            1, static_cast<int>(state.stride.a[0]),
            static_cast<int>(state.stride.a[1])};

        FOR (icomp, jcomp)
        {
            const int ivar      = ivar_0 + VAR_IDX0(icomp, jcomp);
            const auto *var_ptr = state_ptr_xyz + ivar * state.stride.a[2];
            FOR (idir)
            {
                d1(icomp, jcomp, idir) = diff1(var_ptr, strides[idir]);
            }
        }
        return d1;
    }

    // gets the derivative of a consecutive series of vars in a state
    template <int num_diff_vars>
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        amrex::Array2D<amrex::Real, 0, num_diff_vars - 1, 0, AMREX_SPACEDIM - 1>
        diff1_state(int ix, int iy, int iz,
                    const amrex::Array4<const amrex::Real> &state,
                    int first_var = 0) const
    {
        const auto *state_ptr_xyz = state.ptr(ix, iy, iz);
        amrex::GpuArray<int, AMREX_SPACEDIM> strides{
            1, static_cast<int>(state.stride.a[0]),
            static_cast<int>(state.stride.a[1])};

        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
        amrex::Array2D<amrex::Real, 0, num_diff_vars - 1, 0, AMREX_SPACEDIM - 1>
            d1_state{};

        for (int ivar = first_var; ivar < (first_var + num_diff_vars); ivar++)
        {
            const auto *var_ptr = state_ptr_xyz + ivar * state.stride.a[2];

            FOR (idir)
            {
                d1_state(ivar, idir) = diff1(var_ptr, strides[idir]);
            }
        }
        return d1_state;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    diff2(const amrex::Real *in_ptr, const int stride, const int idx = 0) const
    {
        amrex::Real weight_far   = 8.33333333333333333333e-2_rt;
        amrex::Real weight_near  = 1.33333333333333333333e+0_rt;
        amrex::Real weight_local = 2.50000000000000000000e+0_rt;

        return (-weight_far * in_ptr[idx - 2 * stride] +
                weight_near * in_ptr[idx - stride] -
                weight_local * in_ptr[idx] +
                weight_near * in_ptr[idx + stride] -
                weight_far * in_ptr[idx + 2 * stride]) *
               m_one_over_dx2;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    mixed_diff2(const amrex::Real *in_ptr, const int stride1, const int stride2,
                const int idx = 0) const
    {
        amrex::Real weight_far_far   = 6.94444444444444444444e-3_rt;
        amrex::Real weight_near_far  = 5.55555555555555555556e-2_rt;
        amrex::Real weight_near_near = 4.44444444444444444444e-1_rt;

        return (weight_far_far * in_ptr[idx - 2 * stride1 - 2 * stride2] -
                weight_near_far * in_ptr[idx - 2 * stride1 - stride2] +
                weight_near_far * in_ptr[idx - 2 * stride1 + stride2] -
                weight_far_far * in_ptr[idx - 2 * stride1 + 2 * stride2]

                - weight_near_far * in_ptr[idx - stride1 - 2 * stride2] +
                weight_near_near * in_ptr[idx - stride1 - stride2] -
                weight_near_near * in_ptr[idx - stride1 + stride2] +
                weight_near_far * in_ptr[idx - stride1 + 2 * stride2]

                + weight_near_far * in_ptr[idx + stride1 - 2 * stride2] -
                weight_near_near * in_ptr[idx + stride1 - stride2] +
                weight_near_near * in_ptr[idx + stride1 + stride2] -
                weight_near_far * in_ptr[idx + stride1 + 2 * stride2]

                - weight_far_far * in_ptr[idx + 2 * stride1 - 2 * stride2] +
                weight_near_far * in_ptr[idx + 2 * stride1 - stride2] -
                weight_near_far * in_ptr[idx + 2 * stride1 + stride2] +
                weight_far_far * in_ptr[idx + 2 * stride1 + 2 * stride2]) *
               m_one_over_dx2;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE TensorArray::Rank1Sym
    diff2_scalar(int ix, int iy, int iz,
                 const amrex::Array4<amrex::Real const> &state,
                 const int ivar) const
    {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
        TensorArray::Rank1Sym d2;
        const auto *state_ptr_xyz = state.ptr(ix, iy, iz);
        amrex::GpuArray<int, AMREX_SPACEDIM> strides{
            1, static_cast<int>(state.stride.a[0]),
            static_cast<int>(state.stride.a[1])};
        const auto *var_ptr = state_ptr_xyz + ivar * state.stride.a[2];

        d2(0) = diff2(var_ptr, strides[0]);
        d2(3) = diff2(var_ptr, strides[1]);
        d2(5) = diff2(var_ptr, strides[2]);

        d2(1) = mixed_diff2(var_ptr, strides[0], strides[1]);
        d2(2) = mixed_diff2(var_ptr, strides[0], strides[2]);
        d2(4) = mixed_diff2(var_ptr, strides[1], strides[2]);

        return d2;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        amrex::Array2D<amrex::Real, 0, AMREX_SPACEDIM - 1, 0, NUM_SYM_IDXS - 1>
        diff2_vector(int ix, int iy, int iz,
                     const amrex::Array4<amrex::Real const> &state,
                     const int ivar_0) const
    {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
        amrex::Array2D<amrex::Real, 0, AMREX_SPACEDIM - 1, 0, NUM_SYM_IDXS - 1>
            d2;
        const auto *state_ptr_xyz = state.ptr(ix, iy, iz);
        amrex::GpuArray<int, AMREX_SPACEDIM> strides{
            1, static_cast<int>(state.stride.a[0]),
            static_cast<int>(state.stride.a[1])};

        FOR (icomp)
        {
            const int ivar      = ivar_0 + icomp;
            const auto *var_ptr = state_ptr_xyz + ivar * state.stride.a[2];

            d2(icomp, 0) = diff2(var_ptr, strides[0]);
            d2(icomp, 3) = diff2(var_ptr, strides[1]);
            d2(icomp, 5) = diff2(var_ptr, strides[2]);

            d2(icomp, 1) = mixed_diff2(var_ptr, strides[0], strides[1]);
            d2(icomp, 2) = mixed_diff2(var_ptr, strides[0], strides[2]);
            d2(icomp, 4) = mixed_diff2(var_ptr, strides[1], strides[2]);
        }
        return d2;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE TensorArray::Rank2Sym
    diff2_tensor(int ix, int iy, int iz,
                 const amrex::Array4<amrex::Real const> &state,
                 const int ivar_0) const
    {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
        TensorArray::Rank2Sym d2;
        const auto *state_ptr_xyz = state.ptr(ix, iy, iz);
        amrex::GpuArray<int, AMREX_SPACEDIM> strides{
            1, static_cast<int>(state.stride.a[0]),
            static_cast<int>(state.stride.a[1])};

        //        FOR (icomp)
        for (int icomp = 0; icomp < NUM_SYM_IDXS; ++icomp)
        {
            const int ivar      = ivar_0 + icomp;
            const auto *var_ptr = state_ptr_xyz + ivar * state.stride.a[2];

            d2(icomp, 0) = diff2(var_ptr, strides[0]);
            d2(icomp, 3) = diff2(var_ptr, strides[1]);
            d2(icomp, 5) = diff2(var_ptr, strides[2]);

            d2(icomp, 1) = mixed_diff2(var_ptr, strides[0], strides[1]);
            d2(icomp, 2) = mixed_diff2(var_ptr, strides[0], strides[2]);
            d2(icomp, 4) = mixed_diff2(var_ptr, strides[1], strides[2]);
        }

        return d2;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor::Sym12Rank2
    diff2_scalar_test(int ix, int iy, int iz,
                      const amrex::Array4<amrex::Real const> &state,
                      const int ivar) const
    {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
        Tensor::Sym12Rank2 d2;
        const auto *state_ptr_xyz = state.ptr(ix, iy, iz);
        amrex::GpuArray<int, AMREX_SPACEDIM> strides{
            1, static_cast<int>(state.stride.a[0]),
            static_cast<int>(state.stride.a[1])};
        const auto *var_ptr = state_ptr_xyz + ivar * state.stride.a[2];

        d2(0) = diff2(var_ptr, strides[0]);
        d2(3) = diff2(var_ptr, strides[1]);
        d2(5) = diff2(var_ptr, strides[2]);

        d2(1) = mixed_diff2(var_ptr, strides[0], strides[1]);
        d2(2) = mixed_diff2(var_ptr, strides[0], strides[2]);
        d2(4) = mixed_diff2(var_ptr, strides[1], strides[2]);

        return d2;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor::Sym23Rank3
    diff2_vector_test(int ix, int iy, int iz,
                      const amrex::Array4<amrex::Real const> &state,
                      const int ivar_0) const
    {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
        Tensor::Sym23Rank3 d2;
        const auto *state_ptr_xyz = state.ptr(ix, iy, iz);
        amrex::GpuArray<int, AMREX_SPACEDIM> strides{
            1, static_cast<int>(state.stride.a[0]),
            static_cast<int>(state.stride.a[1])};

        FOR (icomp)
        {
            const int ivar      = ivar_0 + icomp;
            const auto *var_ptr = state_ptr_xyz + ivar * state.stride.a[2];

            d2(icomp, 0) = diff2(var_ptr, strides[0]);
            d2(icomp, 3) = diff2(var_ptr, strides[1]);
            d2(icomp, 5) = diff2(var_ptr, strides[2]);

            d2(icomp, 1) = mixed_diff2(var_ptr, strides[0], strides[1]);
            d2(icomp, 2) = mixed_diff2(var_ptr, strides[0], strides[2]);
            d2(icomp, 4) = mixed_diff2(var_ptr, strides[1], strides[2]);
        }
        return d2;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor::Sym12Sym34Rank4
    diff2_tensor_test(int ix, int iy, int iz,
                      const amrex::Array4<amrex::Real const> &state,
                      const int ivar_0) const
    {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
        Tensor::Sym12Sym34Rank4 d2;
        const auto *state_ptr_xyz = state.ptr(ix, iy, iz);
        amrex::GpuArray<int, AMREX_SPACEDIM> strides{
            1, static_cast<int>(state.stride.a[0]),
            static_cast<int>(state.stride.a[1])};

        //        FOR (icomp)
        for (int icomp = 0; icomp < NUM_SYM_IDXS; ++icomp)
        {
            const int ivar      = ivar_0 + icomp;
            const auto *var_ptr = state_ptr_xyz + ivar * state.stride.a[2];

            d2(icomp, 0) = diff2(var_ptr, strides[0]);
            d2(icomp, 3) = diff2(var_ptr, strides[1]);
            d2(icomp, 5) = diff2(var_ptr, strides[2]);

            d2(icomp, 1) = mixed_diff2(var_ptr, strides[0], strides[1]);
            d2(icomp, 2) = mixed_diff2(var_ptr, strides[0], strides[2]);
            d2(icomp, 4) = mixed_diff2(var_ptr, strides[1], strides[2]);
        }

        return d2;
    }

  protected:
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    advection_term(const amrex::Real *in_ptr, const amrex::Real &shift_comp,
                   const int stride, const bool shift_positive,
                   const int idx = 0) const
    {
        amrex::Real weight_0 = -2.50000000000000000000e-1_rt;
        amrex::Real weight_1 = -8.33333333333333333333e-1_rt;
        amrex::Real weight_2 = +1.50000000000000000000e+0_rt;
        amrex::Real weight_3 = -5.00000000000000000000e-1_rt;
        amrex::Real weight_4 = +8.33333333333333333333e-2_rt;

        amrex::Real upwind =
            shift_comp *
            (weight_0 * in_ptr[idx - stride] + weight_1 * in_ptr[idx] +
             weight_2 * in_ptr[idx + stride] +
             weight_3 * in_ptr[idx + 2 * stride] +
             weight_4 * in_ptr[idx + 3 * stride]) *
            m_one_over_dx;

        amrex::Real downwind =
            shift_comp *
            (-weight_4 * in_ptr[idx - 3 * stride] -
             weight_3 * in_ptr[idx - 2 * stride] -
             weight_2 * in_ptr[idx - stride] - weight_1 * in_ptr[idx] -
             weight_0 * in_ptr[idx + stride]) *
            m_one_over_dx;

        return (shift_positive) ? upwind : downwind;
    }

  public:

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    advection(int ix, int iy, int iz,
              const amrex::Array4<amrex::Real const> &state,
              const TensorArray::Rank1 &shift_vector, const int ivar) const
    {
        amrex::Real advec         = 0.0;
        const auto *state_ptr_xyz = state.ptr(ix, iy, iz);
        int j_stride              = static_cast<int>(state.stride.a[0]);
        int k_stride              = static_cast<int>(state.stride.a[1]);
        int n_stride              = static_cast<int>(state.stride.a[2]);

        amrex::GpuArray<int, AMREX_SPACEDIM> strides{1, j_stride, k_stride};

        const auto *var_ptr =
            state_ptr_xyz + static_cast<amrex::Long>(ivar) * n_stride;
        FOR (idir)
        {
            const bool shift_positive = (shift_vector(idir) > 0.0);
            advec += advection_term(var_ptr, shift_vector(idir), strides[idir],
                                    shift_positive);
        }
        return advec;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    advection(int ix, int iy, int iz,
              const amrex::Array4<amrex::Real const> &state,
              const Tensor::Rank1 &shift_vector, const int ivar) const
    {
        amrex::Real advec         = 0.0;
        const auto *state_ptr_xyz = state.ptr(ix, iy, iz);
        int j_stride              = static_cast<int>(state.stride.a[0]);
        int k_stride              = static_cast<int>(state.stride.a[1]);
        int n_stride              = static_cast<int>(state.stride.a[2]);

        amrex::GpuArray<int, AMREX_SPACEDIM> strides{1, j_stride, k_stride};

        const auto *var_ptr =
            state_ptr_xyz + static_cast<amrex::Long>(ivar) * n_stride;
        FOR (idir)
        {
            const bool shift_positive = (shift_vector(idir) > 0.0);
            advec += advection_term(var_ptr, shift_vector(idir), strides[idir],
                                    shift_positive);
        }
        return advec;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    advec_scalar(int ix, int iy, int iz,
                 const amrex::Array4<amrex::Real const> &state,
                 const TensorArray::Rank1 &shift_vector, const int ivar) const
    {
        return advection(ix, iy, iz, state, shift_vector, ivar);
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    advec_scalar(int ix, int iy, int iz,
                 const amrex::Array4<amrex::Real const> &state,
                 const Tensor::Rank1 &shift_vector, const int ivar) const
    {
        return advection(ix, iy, iz, state, shift_vector, ivar);
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE TensorArray::Rank1
    advec_vector(int ix, int iy, int iz,
                 const amrex::Array4<amrex::Real const> &state,
                 const TensorArray::Rank1 &shift_vector, const int ivar0) const
    {
        TensorArray::Rank1 advec_vector{};
        FOR (icomp)
        {
            int ivar = ivar0 + icomp;
            advec_vector(icomp) =
                advection(ix, iy, iz, state, shift_vector, ivar);
        }
        return advec_vector;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor::Rank1
    advec_vector(int ix, int iy, int iz,
                 const amrex::Array4<amrex::Real const> &state,
                 const Tensor::Rank1 &shift_vector, const int ivar0) const
    {
        Tensor::Rank1 advec_vector{};
        FOR (icomp)
        {
            int ivar = ivar0 + icomp;
            advec_vector(icomp) =
                advection(ix, iy, iz, state, shift_vector, ivar);
        }
        return advec_vector;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE TensorArray::Rank2
    advec_tensor(int ix, int iy, int iz,
                 const amrex::Array4<amrex::Real const> &state,
                 const TensorArray::Rank1 &shift_vector, const int ivar0) const
    {
        TensorArray::Rank2 advec_tensor{};
        FOR (icomp, jcomp)
        {
            int ivar = VAR_IDX(ivar0, icomp, jcomp);
            advec_tensor(icomp, jcomp) =
                advection(ix, iy, iz, state, shift_vector, ivar);
        }
        return advec_tensor;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor::Rank2
    advec_tensor(int ix, int iy, int iz,
                 const amrex::Array4<amrex::Real const> &state,
                 const Tensor::Rank1 &shift_vector, const int ivar0) const
    {
        Tensor::Rank2 advec_tensor{};
        FOR (icomp, jcomp)
        {
            int ivar = VAR_IDX(ivar0, icomp, jcomp);
            advec_tensor(icomp, jcomp) =
                advection(ix, iy, iz, state, shift_vector, ivar);
        }
        return advec_tensor;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE TensorArray::Rank1Sym
    advec_sym_tensor(int ix, int iy, int iz,
                     const amrex::Array4<amrex::Real const> &state,
                     const TensorArray::Rank1 &shift_vector,
                     const int ivar0) const
    {
        TensorArray::Rank1Sym advec_tensor{};

        for (int i = 0; i < TensorArray::Rank1Sym::len(); i++)
        {
            advec_tensor(i) =
                advection(ix, iy, iz, state, shift_vector, ivar0 + i);
        }
        return advec_tensor;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor::Sym12Rank2
    advec_sym_tensor(int ix, int iy, int iz,
                     const amrex::Array4<amrex::Real const> &state,
                     const Tensor::Rank1 &shift_vector, const int ivar0) const
    {
        Tensor::Sym12Rank2 advec_tensor{};

        for (int i = 0; i < NUM_SYM_IDXS; i++)
        {
            advec_tensor(i) =
                advection(ix, iy, iz, state, shift_vector, ivar0 + i);
        }
        return advec_tensor;
    }

    // gets the derivative of a consecutive series of vars in a state
    template <int num_diff_vars>
    [[nodiscard]] AMREX_GPU_DEVICE
        AMREX_FORCE_INLINE amrex::GpuArray<amrex::Real, num_diff_vars>
        advec_state(int ix, int iy, int iz,
                    const amrex::Array4<const amrex::Real> &state,
                    const TensorArray::Rank1 &shift_vector,
                    int first_var = 0) const

    {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
        amrex::GpuArray<amrex::Real, num_diff_vars> advec_state;
        for (int ivar = first_var; ivar < (first_var + num_diff_vars); ivar++)
        {
            advec_state[ivar] =
                advection(ix, iy, iz, state, shift_vector, ivar);
        }
        return advec_state;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    dissipation_term(const double *in_ptr, const int stride,
                     const int idx = 0) const
    {
        amrex::Real weight_vfar  = 1.56250e-2_rt;
        amrex::Real weight_far   = 9.37500e-2_rt;
        amrex::Real weight_near  = 2.34375e-1_rt;
        amrex::Real weight_local = 3.12500e-1_rt;

        return (weight_vfar * in_ptr[idx - 3 * stride] -
                weight_far * in_ptr[idx - 2 * stride] +
                weight_near * in_ptr[idx - stride] -
                weight_local * in_ptr[idx] +
                weight_near * in_ptr[idx + stride] -
                weight_far * in_ptr[idx + 2 * stride] +
                weight_vfar * in_ptr[idx + 3 * stride]) *
               m_one_over_dx;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    calculate_dissipation(int ix, int iy, int iz,
                          const amrex::Array4<amrex::Real const> &state,
                          const double sigma_coeff, const int ivar) const
    {
        amrex::Real diss          = 0.0;
        const auto *state_ptr_xyz = state.ptr(ix, iy, iz);

        int j_stride = static_cast<int>(state.stride.a[0]);
        int k_stride = static_cast<int>(state.stride.a[1]);
        int n_stride = static_cast<int>(state.stride.a[2]);

        amrex::GpuArray<int, AMREX_SPACEDIM> strides{1, j_stride, k_stride};

        FOR (idir)
        {
            const auto stride = strides[idir];
            diss +=
                sigma_coeff *
                dissipation_term(state_ptr_xyz +
                                     static_cast<amrex::Long>(ivar) * n_stride,
                                 stride);
        }
        return diss;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    add_dissipation(int ix, int iy, int iz,
                    const amrex::CellData<amrex::Real> &rhs_cell_data,
                    const amrex::Array4<amrex::Real const> &state,
                    const double sigma_coeff, int num_vars = NUM_VARS) const
    {
        for (int ivar = 0; ivar < num_vars; ivar++)
        {
            amrex::Real diss =
                calculate_dissipation(ix, iy, iz, state, sigma_coeff, ivar);
            rhs_cell_data[ivar] += diss;
        }
    }

    // NOLINTEND(bugprone-easily-swappable-parameters)
};

#endif /* FOURTHORDERDERIVATIVES_HPP_ */
