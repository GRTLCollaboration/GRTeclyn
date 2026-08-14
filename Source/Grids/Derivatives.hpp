/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef DERIVATIVES_HPP_
#define DERIVATIVES_HPP_

#include "DimensionDefinitions.hpp"

#include "StateVariables.hpp"
#include "Tensor.hpp"
#include <AMReX_REAL.H>
#include <array>
#include <concepts>
#include "AMReX_Array.H"

using namespace amrex::literals;

/// Concept for the order-specific stencil types used by the Derivatives
/// class. A stencil_t provides the finite difference weights for the pointwise
/// derivative operators (first, second, mixed second, advection and
/// dissipation).
///
/// Note: the signatures must match exactly (the const member function/noexcept
/// qualifier is part of the signature).
template <class stencil_t>
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
concept Stencil =
    requires(const stencil_t stencil, const amrex::Real *in_ptr,
             const amrex::Real &shift_comp, const int stride, const int stride2,
             const bool shift_positive, const int idx)
{
    {
        stencil.diff1(in_ptr, stride, idx)
    } noexcept -> std::same_as<amrex::Real>;
    {
        stencil.diff2(in_ptr, stride, idx)
    } noexcept -> std::same_as<amrex::Real>;
    {
        stencil.mixed_diff2(in_ptr, stride, stride2, idx)
    } noexcept -> std::same_as<amrex::Real>;
    {
        stencil.advection_term(in_ptr, shift_comp, stride, shift_positive, idx)
    } noexcept -> std::same_as<amrex::Real>;
    {
        stencil.dissipation_term(in_ptr, stride, idx)
    } noexcept -> std::same_as<amrex::Real>;
    // The sign applied to the dissipation term in calculate_dissipation
    {
        stencil_t::dissipation_sign
    } noexcept -> std::convertible_to<amrex::Real>;
};

// NOLINTBEGIN(bugprone-easily-swappable-parameters)

/// Fourth order stencils (with sixth order KO dissipation)
struct FourthOrderStencil
{
    /// Sign of the dissipation in calculate_dissipation (plus for fourth
    /// order/sixth order dissipation, minus for sixth order/eighth order
    /// dissipation)
    static constexpr amrex::Real dissipation_sign = 1.0;

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE static amrex::Real
    diff1(const amrex::Real *in_ptr, const int stride,
          const int idx = 0) noexcept
    {
        amrex::Real weight_far  = 8.33333333333333333333e-2_rt;
        amrex::Real weight_near = 6.66666666666666666667e-1_rt;

        return weight_far * in_ptr[idx - 2 * stride] -
               weight_near * in_ptr[idx - stride] +
               weight_near * in_ptr[idx + stride] -
               weight_far * in_ptr[idx + 2 * stride];
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE static amrex::Real
    diff2(const amrex::Real *in_ptr, const int stride,
          const int idx = 0) noexcept
    {
        amrex::Real weight_far   = 8.33333333333333333333e-2_rt;
        amrex::Real weight_near  = 1.33333333333333333333e+0_rt;
        amrex::Real weight_local = 2.50000000000000000000e+0_rt;

        return -weight_far * in_ptr[idx - 2 * stride] +
               weight_near * in_ptr[idx - stride] - weight_local * in_ptr[idx] +
               weight_near * in_ptr[idx + stride] -
               weight_far * in_ptr[idx + 2 * stride];
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE static amrex::Real
    mixed_diff2(const amrex::Real *in_ptr, const int stride1, const int stride2,
                const int idx = 0) noexcept
    {
        amrex::Real weight_far_far   = 6.94444444444444444444e-3_rt;
        amrex::Real weight_near_far  = 5.55555555555555555556e-2_rt;
        amrex::Real weight_near_near = 4.44444444444444444444e-1_rt;

        return weight_far_far * in_ptr[idx - 2 * stride1 - 2 * stride2] -
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
               weight_far_far * in_ptr[idx + 2 * stride1 + 2 * stride2];
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE static amrex::Real
    advection_term(const amrex::Real *in_ptr, const amrex::Real &shift_comp,
                   const int stride, const bool shift_positive,
                   const int idx = 0) noexcept
    {
        amrex::Real weight_0 = -2.50000000000000000000e-1_rt;
        amrex::Real weight_1 = -8.33333333333333333333e-1_rt;
        amrex::Real weight_2 = +1.50000000000000000000e+0_rt;
        amrex::Real weight_3 = -5.00000000000000000000e-1_rt;
        amrex::Real weight_4 = +8.33333333333333333333e-2_rt;

        amrex::Real upwind = shift_comp * (weight_0 * in_ptr[idx - stride] +
                                           weight_1 * in_ptr[idx] +
                                           weight_2 * in_ptr[idx + stride] +
                                           weight_3 * in_ptr[idx + 2 * stride] +
                                           weight_4 * in_ptr[idx + 3 * stride]);

        amrex::Real downwind =
            shift_comp *
            (-weight_4 * in_ptr[idx - 3 * stride] -
             weight_3 * in_ptr[idx - 2 * stride] -
             weight_2 * in_ptr[idx - stride] - weight_1 * in_ptr[idx] -
             weight_0 * in_ptr[idx + stride]);

        return (shift_positive) ? upwind : downwind;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE static amrex::Real
    dissipation_term(const amrex::Real *in_ptr, const int stride,
                     const int idx = 0) noexcept
    {
        amrex::Real weight_vfar  = 1.56250e-2_rt;
        amrex::Real weight_far   = 9.37500e-2_rt;
        amrex::Real weight_near  = 2.34375e-1_rt;
        amrex::Real weight_local = 3.12500e-1_rt;

        return weight_vfar * in_ptr[idx - 3 * stride] -
               weight_far * in_ptr[idx - 2 * stride] +
               weight_near * in_ptr[idx - stride] - weight_local * in_ptr[idx] +
               weight_near * in_ptr[idx + stride] -
               weight_far * in_ptr[idx + 2 * stride] +
               weight_vfar * in_ptr[idx + 3 * stride];
    }
};

/// Sixth order stencils (with eighth order KO dissipation)
struct SixthOrderStencil
{
    /// Sign of the dissipation in calculate_dissipation (minus for sixth
    /// order/eighth order dissipation)
    static constexpr amrex::Real dissipation_sign = -1.0;

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE static amrex::Real
    diff1(const amrex::Real *in_ptr, const int stride,
          const int idx = 0) noexcept
    {
        amrex::Real weight_vfar = 1.66666666666666666667e-2_rt;
        amrex::Real weight_far  = 1.50000000000000000000e-1_rt;
        amrex::Real weight_near = 7.50000000000000000000e-1_rt;

        return -weight_vfar * in_ptr[idx - 3 * stride] +
               weight_far * in_ptr[idx - 2 * stride] -
               weight_near * in_ptr[idx - stride] +
               weight_near * in_ptr[idx + stride] -
               weight_far * in_ptr[idx + 2 * stride] +
               weight_vfar * in_ptr[idx + 3 * stride];
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE static amrex::Real
    diff2(const amrex::Real *in_ptr, const int stride,
          const int idx = 0) noexcept
    {
        amrex::Real weight_vfar  = 1.11111111111111111111e-2_rt;
        amrex::Real weight_far   = 1.50000000000000000000e-1_rt;
        amrex::Real weight_near  = 1.50000000000000000000e+0_rt;
        amrex::Real weight_local = 2.72222222222222222222e+0_rt;

        return weight_vfar * in_ptr[idx - 3 * stride] -
               weight_far * in_ptr[idx - 2 * stride] +
               weight_near * in_ptr[idx - stride] - weight_local * in_ptr[idx] +
               weight_near * in_ptr[idx + stride] -
               weight_far * in_ptr[idx + 2 * stride] +
               weight_vfar * in_ptr[idx + 3 * stride];
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE static amrex::Real
    mixed_diff2(const amrex::Real *in_ptr, const int stride1, const int stride2,
                const int idx = 0) noexcept
    {
        amrex::Real weight_vfar_vfar = 2.77777777777777777778e-4_rt;
        amrex::Real weight_vfar_far  = 2.50000000000000000000e-3_rt;
        amrex::Real weight_vfar_near = 1.25000000000000000000e-2_rt;
        amrex::Real weight_far_far   = 2.25000000000000000000e-2_rt;
        amrex::Real weight_far_near  = 1.12500000000000000000e-1_rt;
        amrex::Real weight_near_near = 5.62500000000000000000e-1_rt;

        return weight_vfar_vfar * in_ptr[idx - 3 * stride1 - 3 * stride2] -
               weight_vfar_far * in_ptr[idx - 3 * stride1 - 2 * stride2] +
               weight_vfar_near * in_ptr[idx - 3 * stride1 - stride2] -
               weight_vfar_near * in_ptr[idx - 3 * stride1 + stride2] +
               weight_vfar_far * in_ptr[idx - 3 * stride1 + 2 * stride2] -
               weight_vfar_vfar * in_ptr[idx - 3 * stride1 + 3 * stride2]

               - weight_vfar_far * in_ptr[idx - 2 * stride1 - 3 * stride2] +
               weight_far_far * in_ptr[idx - 2 * stride1 - 2 * stride2] -
               weight_far_near * in_ptr[idx - 2 * stride1 - stride2] +
               weight_far_near * in_ptr[idx - 2 * stride1 + stride2] -
               weight_far_far * in_ptr[idx - 2 * stride1 + 2 * stride2] +
               weight_vfar_far * in_ptr[idx - 2 * stride1 + 3 * stride2]

               + weight_vfar_near * in_ptr[idx - stride1 - 3 * stride2] -
               weight_far_near * in_ptr[idx - stride1 - 2 * stride2] +
               weight_near_near * in_ptr[idx - stride1 - stride2] -
               weight_near_near * in_ptr[idx - stride1 + stride2] +
               weight_far_near * in_ptr[idx - stride1 + 2 * stride2] -
               weight_vfar_near * in_ptr[idx - stride1 + 3 * stride2]

               - weight_vfar_near * in_ptr[idx + stride1 - 3 * stride2] +
               weight_far_near * in_ptr[idx + stride1 - 2 * stride2] -
               weight_near_near * in_ptr[idx + stride1 - stride2] +
               weight_near_near * in_ptr[idx + stride1 + stride2] -
               weight_far_near * in_ptr[idx + stride1 + 2 * stride2] +
               weight_vfar_near * in_ptr[idx + stride1 + 3 * stride2]

               + weight_vfar_far * in_ptr[idx + 2 * stride1 - 3 * stride2] -
               weight_far_far * in_ptr[idx + 2 * stride1 - 2 * stride2] +
               weight_far_near * in_ptr[idx + 2 * stride1 - stride2] -
               weight_far_near * in_ptr[idx + 2 * stride1 + stride2] +
               weight_far_far * in_ptr[idx + 2 * stride1 + 2 * stride2] -
               weight_vfar_far * in_ptr[idx + 2 * stride1 + 3 * stride2]

               - weight_vfar_vfar * in_ptr[idx + 3 * stride1 - 3 * stride2] +
               weight_vfar_far * in_ptr[idx + 3 * stride1 - 2 * stride2] -
               weight_vfar_near * in_ptr[idx + 3 * stride1 - stride2] +
               weight_vfar_near * in_ptr[idx + 3 * stride1 + stride2] -
               weight_vfar_far * in_ptr[idx + 3 * stride1 + 2 * stride2] +
               weight_vfar_vfar * in_ptr[idx + 3 * stride1 + 3 * stride2];
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE static amrex::Real
    advection_term(const amrex::Real *in_ptr, const amrex::Real &shift_comp,
                   const int stride, const bool shift_positive,
                   const int idx = 0) noexcept
    {
        amrex::Real weight_0 = +3.33333333333333333333e-2_rt;
        amrex::Real weight_1 = -4.00000000000000000000e-1_rt;
        amrex::Real weight_2 = -5.83333333333333333333e-1_rt;
        amrex::Real weight_3 = +1.33333333333333333333e+0_rt;
        amrex::Real weight_4 = -5.00000000000000000000e-1_rt;
        amrex::Real weight_5 = +1.33333333333333333333e-1_rt;
        amrex::Real weight_6 = -1.66666666666666666667e-2_rt;

        amrex::Real upwind = shift_comp * (weight_0 * in_ptr[idx - 2 * stride] +
                                           weight_1 * in_ptr[idx - stride] +
                                           weight_2 * in_ptr[idx] +
                                           weight_3 * in_ptr[idx + stride] +
                                           weight_4 * in_ptr[idx + 2 * stride] +
                                           weight_5 * in_ptr[idx + 3 * stride] +
                                           weight_6 * in_ptr[idx + 4 * stride]);

        amrex::Real downwind =
            shift_comp *
            (-weight_6 * in_ptr[idx - 4 * stride] -
             weight_5 * in_ptr[idx - 3 * stride] -
             weight_4 * in_ptr[idx - 2 * stride] -
             weight_3 * in_ptr[idx - stride] - weight_2 * in_ptr[idx] -
             weight_1 * in_ptr[idx + stride] -
             weight_0 * in_ptr[idx + 2 * stride]);

        return (shift_positive) ? upwind : downwind;
    }

    // Eighth order dissipation, note that this changes the sign of
    // dissipation_sign relative to the fourth order stencil
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE static amrex::Real
    dissipation_term(const amrex::Real *in_ptr, const int stride,
                     const int idx = 0) noexcept
    {
        amrex::Real weight_vvfar = 3.906250e-3_rt;
        amrex::Real weight_vfar  = 3.125000e-2_rt;
        amrex::Real weight_far   = 1.093750e-1_rt;
        amrex::Real weight_near  = 2.187500e-1_rt;
        amrex::Real weight_local = 2.734375e-1_rt;

        return weight_vvfar * in_ptr[idx - 4 * stride] -
               weight_vfar * in_ptr[idx - 3 * stride] +
               weight_far * in_ptr[idx - 2 * stride] -
               weight_near * in_ptr[idx - stride] + weight_local * in_ptr[idx] -
               weight_near * in_ptr[idx + stride] +
               weight_far * in_ptr[idx + 2 * stride] -
               weight_vfar * in_ptr[idx + 3 * stride] +
               weight_vvfar * in_ptr[idx + 4 * stride];
    }
};

static_assert(Stencil<FourthOrderStencil>);
static_assert(Stencil<SixthOrderStencil>);

/// The general derivatives class conforming to the deriv_t API. All of the
/// order-specific finite difference weights are provided by the stencil_t
/// class (see the Stencil concept above).
template <Stencil stencil_t> class Derivatives
{
  protected:
    amrex::Real m_dx;
    amrex::Real m_one_over_dx;
    amrex::Real m_one_over_dx2;

  public:
    AMREX_GPU_HOST_DEVICE
    Derivatives(amrex::Real dx)
        : m_dx(dx), m_one_over_dx(1 / dx), m_one_over_dx2(1 / (dx * dx))
    {
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    diff1(const amrex::Real *in_ptr, const int stride, const int idx = 0) const
    {
        return stencil_t::diff1(in_ptr, stride, idx) * m_one_over_dx;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    diff2(const amrex::Real *in_ptr, const int stride, const int idx = 0) const
    {
        return stencil_t::diff2(in_ptr, stride, idx) * m_one_over_dx2;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    mixed_diff2(const amrex::Real *in_ptr, const int stride1, const int stride2,
                const int idx = 0) const
    {
        return stencil_t::mixed_diff2(in_ptr, stride1, stride2, idx) *
               m_one_over_dx2;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    dissipation_term(const amrex::Real *in_ptr, const int stride,
                     const int idx = 0) const
    {
        return stencil_t::dissipation_term(in_ptr, stride, idx) * m_one_over_dx;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor::Rank1
    d1_scalar(int ix, int iy, int iz,
              const amrex::Array4<const amrex::Real> &state,
              const int ivar) const
    {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
        Tensor::Rank1 d1;
        const auto *state_ptr_xyz = state.ptr(ix, iy, iz);
        const auto strides        = get_strides(state);
        const auto *var_ptr       = get_var_ptr(ivar, state_ptr_xyz, strides);
        FOR (idir)
        {
            d1(idir) = diff1(var_ptr, strides[idir]);
        }
        return d1;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor::Rank2
    d1_vector(int ix, int iy, int iz,
              const amrex::Array4<const amrex::Real> &state,
              const int ivar_0) const
    {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
        Tensor::Rank2 d1;
        const auto *state_ptr_xyz = state.ptr(ix, iy, iz);
        const auto strides        = get_strides(state);

        FOR (icomp)
        {
            const int ivar      = ivar_0 + icomp;
            const auto *var_ptr = get_var_ptr(ivar, state_ptr_xyz, strides);
            FOR (idir)
            {
                d1(icomp, idir) = diff1(var_ptr, strides[idir]);
            }
        }
        return d1;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor::Sym12Rank3
    d1_sym_tensor(int ix, int iy, int iz,
                  const amrex::Array4<const amrex::Real> &state,
                  const int ivar_0) const
    {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
        Tensor::Sym12Rank3 d1{};
        const auto *state_ptr_xyz = state.ptr(ix, iy, iz);
        const auto strides        = get_strides(state);

        for (int ivar = 0; ivar < NUM_SYM_IDXS; ++ivar)
        {
            const auto *var_ptr =
                get_var_ptr(ivar_0 + ivar, state_ptr_xyz, strides);

            FOR (idir)
            {
                d1(ivar, idir) = diff1(var_ptr, strides[idir]);
            }
        }
        return d1;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor::Rank3
    d1_tensor(int ix, int iy, int iz,
              const amrex::Array4<const amrex::Real> &state,
              const int ivar_0) const
    {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
        Tensor::Rank3 d1{};
        const auto *state_ptr_xyz = state.ptr(ix, iy, iz);
        const auto strides        = get_strides(state);

        int ivar{ivar_0};
        FOR (icomp, jcomp)
        {
            const auto *var_ptr = get_var_ptr(ivar, state_ptr_xyz, strides);
            FOR (idir)
            {
                d1(icomp, jcomp, idir) = diff1(var_ptr, strides[idir]);
            }
            ++ivar;
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

        for (int ivar = first_var; ivar < (first_var + num_diff_vars); ++ivar)
        {
            const auto *var_ptr = state_ptr_xyz + ivar * state.stride.a[2];

            FOR (idir)
            {
                d1_state(ivar - first_var, idir) =
                    diff1(var_ptr, strides[idir]);
            }
        }
        return d1_state;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor::Sym12Rank2
    d2_scalar(int ix, int iy, int iz,
              const amrex::Array4<const amrex::Real> &state,
              const int ivar) const
    {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
        Tensor::Sym12Rank2 d2;

        const auto *state_ptr_xyz = state.ptr(ix, iy, iz);
        auto strides              = get_strides(state);
        const auto *var_ptr       = get_var_ptr(ivar, state_ptr_xyz, strides);

        d2(0, 0) = diff2(var_ptr, strides[0]);
        d2(1, 1) = diff2(var_ptr, strides[1]);
        d2(2, 2) = diff2(var_ptr, strides[2]);

        d2(0, 1) = mixed_diff2(var_ptr, strides[0], strides[1]);
        d2(0, 2) = mixed_diff2(var_ptr, strides[0], strides[2]);
        d2(1, 2) = mixed_diff2(var_ptr, strides[1], strides[2]);

        return d2;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor::Sym23Rank3
    d2_vector(int ix, int iy, int iz,
              const amrex::Array4<const amrex::Real> &state,
              const int ivar_0) const
    {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
        Tensor::Sym23Rank3 d2;
        const auto *state_ptr_xyz = state.ptr(ix, iy, iz);
        auto strides              = get_strides(state);

        FOR (icomp)
        {
            const int ivar      = ivar_0 + icomp;
            const auto *var_ptr = get_var_ptr(ivar, state_ptr_xyz, strides);

            d2(icomp, 0, 0) = diff2(var_ptr, strides[0]);
            d2(icomp, 1, 1) = diff2(var_ptr, strides[1]);
            d2(icomp, 2, 2) = diff2(var_ptr, strides[2]);

            d2(icomp, 0, 1) = mixed_diff2(var_ptr, strides[0], strides[1]);
            d2(icomp, 0, 2) = mixed_diff2(var_ptr, strides[0], strides[2]);
            d2(icomp, 1, 2) = mixed_diff2(var_ptr, strides[1], strides[2]);
        }
        return d2;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor::Sym34Rank4
    d2_tensor(int ix, int iy, int iz,
              const amrex::Array4<const amrex::Real> &state,
              const int ivar_0) const
    {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
        Tensor::Sym34Rank4 d2;
        const auto *state_ptr_xyz = state.ptr(ix, iy, iz);
        auto strides              = get_strides(state);

        int ivar{ivar_0};
        FOR (icomp)
            FOR (jcomp)
            {
                const auto *var_ptr = get_var_ptr(ivar, state_ptr_xyz, strides);

                d2(icomp, jcomp, 0, 0) = diff2(var_ptr, strides[0]);
                d2(icomp, jcomp, 1, 1) = diff2(var_ptr, strides[1]);
                d2(icomp, jcomp, 2, 2) = diff2(var_ptr, strides[2]);

                d2(icomp, jcomp, 0, 1) =
                    mixed_diff2(var_ptr, strides[0], strides[1]);
                d2(icomp, jcomp, 0, 2) =
                    mixed_diff2(var_ptr, strides[0], strides[2]);
                d2(icomp, jcomp, 1, 2) =
                    mixed_diff2(var_ptr, strides[1], strides[2]);

                ++ivar;
            }

        return d2;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor::Sym12Sym34Rank4
    d2_sym_tensor(int ix, int iy, int iz,
                  const amrex::Array4<const amrex::Real> &state,
                  const int ivar_0) const
    {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
        Tensor::Sym12Sym34Rank4 d2;
        const auto *state_ptr_xyz = state.ptr(ix, iy, iz);
        auto strides              = get_strides(state);

        FOR (icomp)
            FOR (jcomp)
            {
                const int ivar      = sym_var_idx(ivar_0, icomp, jcomp);
                const auto *var_ptr = get_var_ptr(ivar, state_ptr_xyz, strides);

                d2(icomp, jcomp, 0, 0) = diff2(var_ptr, strides[0]);
                d2(icomp, jcomp, 1, 1) = diff2(var_ptr, strides[1]);
                d2(icomp, jcomp, 2, 2) = diff2(var_ptr, strides[2]);

                d2(icomp, jcomp, 0, 1) =
                    mixed_diff2(var_ptr, strides[0], strides[1]);
                d2(icomp, jcomp, 0, 2) =
                    mixed_diff2(var_ptr, strides[0], strides[2]);
                d2(icomp, jcomp, 1, 2) =
                    mixed_diff2(var_ptr, strides[1], strides[2]);
            }

        return d2;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    advection(int ix, int iy, int iz,
              const amrex::Array4<const amrex::Real> &state,
              const Tensor::Rank1 &shift_vector, const int ivar) const
    {
        amrex::Real advec         = 0.0;
        const auto *state_ptr_xyz = state.ptr(ix, iy, iz);
        const auto strides        = get_strides(state);
        const auto *var_ptr       = get_var_ptr(ivar, state_ptr_xyz, strides);

        FOR (idir)
        {
            const bool shift_positive = (shift_vector(idir) > 0.0);
            advec += stencil_t::advection_term(var_ptr, shift_vector(idir),
                                               strides[idir], shift_positive) *
                     m_one_over_dx;
        }
        return advec;
    }

    // NOLINTBEGIN(readability-convert-member-functions-to-static)
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    advec_scalar(int ix, int iy, int iz,
                 const amrex::Array4<const amrex::Real> &state,
                 const Tensor::Rank1 &shift_vector, const int ivar) const
    {
        return advection(ix, iy, iz, state, shift_vector, ivar);
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor::Rank1
    advec_vector(int ix, int iy, int iz,
                 const amrex::Array4<const amrex::Real> &state,
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

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor::Rank2
    advec_tensor(int ix, int iy, int iz,
                 const amrex::Array4<const amrex::Real> &state,
                 const Tensor::Rank1 &shift_vector, const int ivar0) const
    {
        Tensor::Rank2 advec_tensor{};
        FOR (icomp, jcomp)
        {
            int ivar = sym_var_idx(ivar0, icomp, jcomp);
            advec_tensor(icomp, jcomp) =
                advection(ix, iy, iz, state, shift_vector, ivar);
        }
        return advec_tensor;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor::Sym12Rank2
    advec_sym_tensor(int ix, int iy, int iz,
                     const amrex::Array4<const amrex::Real> &state,
                     const Tensor::Rank1 &shift_vector, const int ivar0) const
    {
        Tensor::Sym12Rank2 advec_tensor{};

        for (int i = 0; i < NUM_SYM_IDXS; ++i)
        {
            advec_tensor(i) =
                advection(ix, iy, iz, state, shift_vector, ivar0 + i);
        }
        return advec_tensor;
    }
    // NOLINTEND(readability-convert-member-functions-to-static)

    // gets the derivative of a consecutive series of vars in a state
    template <int num_diff_vars>
    [[nodiscard]] AMREX_GPU_DEVICE
        AMREX_FORCE_INLINE amrex::GpuArray<amrex::Real, num_diff_vars>
        advec_state(int ix, int iy, int iz,
                    const amrex::Array4<const amrex::Real> &state,
                    const Tensor::Rank1 &shift_vector, int first_var = 0) const

    {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
        amrex::GpuArray<amrex::Real, num_diff_vars> advec_state;
        for (int ivar = first_var; ivar < (first_var + num_diff_vars); ++ivar)
        {
            advec_state[ivar - first_var] =
                advection(ix, iy, iz, state, shift_vector, ivar);
        }
        return advec_state;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    calculate_dissipation(int ix, int iy, int iz,
                          const amrex::Array4<const amrex::Real> &state,
                          const amrex::Real sigma_coeff, const int ivar) const
    {
        amrex::Real diss          = 0.0;
        const auto *state_ptr_xyz = state.ptr(ix, iy, iz);
        const auto strides        = get_strides(state);

        FOR (idir)
        {
            const auto stride  = strides[idir];
            diss              += stencil_t::dissipation_sign * sigma_coeff *
                    dissipation_term(get_var_ptr(ivar, state_ptr_xyz, strides),
                                     stride);
        }
        return diss;
    }

    // NOLINTBEGIN(readability-convert-member-functions-to-static)
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    add_dissipation(int ix, int iy, int iz,
                    const amrex::CellData<amrex::Real> &rhs_cell_data,
                    const amrex::Array4<const amrex::Real> &state,
                    const amrex::Real sigma_coeff,
                    int num_vars = NUM_VARS) const

    {
        amrex::Real diss = 0.0;
        for (int ivar = 0; ivar < num_vars; ++ivar)
        {
            diss = calculate_dissipation(ix, iy, iz, state, sigma_coeff, ivar);
            rhs_cell_data[ivar] += diss;
        }
    }
    // NOLINTEND(readability-convert-member-functions-to-static)

  protected:
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

// NOLINTEND(bugprone-easily-swappable-parameters)

using FourthOrderDerivatives = Derivatives<FourthOrderStencil>;
using SixthOrderDerivatives  = Derivatives<SixthOrderStencil>;

#endif /* DERIVATIVES_HPP_ */
