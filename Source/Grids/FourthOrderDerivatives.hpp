/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef FOURTHORDERDERIVATIVES_HPP_
#define FOURTHORDERDERIVATIVES_HPP_

#include "Cell.hpp"
#include "DimensionDefinitions.hpp"
#include <AMReX_Array.H>
#include <AMReX_REAL.H>
#include <array>

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

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    diff1(const amrex::Real *in_ptr, const int idx, const int stride) const
    {
        amrex::Real weight_far  = 8.33333333333333333333e-2_rt;
        amrex::Real weight_near = 6.66666666666666666667e-1_rt;

        return (weight_far * in_ptr[idx - 2 * stride] -
                weight_near * in_ptr[idx - stride] +
                weight_near * in_ptr[idx + stride] -
                weight_far * in_ptr[idx + 2 * stride]) *
               m_one_over_dx;
    }

    /// Calculates all first derivatives and returns as variable type specified
    /// by the template parameter
    template <template <typename> class vars_t>
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE auto
    diff1(int i, int j, int k,
          const amrex::Array4<const amrex::Real> &state) const
    {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
        vars_t<amrex::Array1D<amrex::Real, 0, AMREX_SPACEDIM>> d1;
        const auto *state_ptr_ijk = state.ptr(i, j, k);
        d1.enum_mapping(
            [&](const int &ivar,
                amrex::Array1D<amrex::Real, 0, AMREX_SPACEDIM> &var)
            {
                AMREX_D_TERM(
                    var(0) = diff1(state_ptr_ijk + ivar * state.nstride, 0, 1);
                    , var(1) = diff1(state_ptr_ijk + ivar * state.nstride, 0,
                                     static_cast<int>(state.jstride));
                    , var(2) = diff1(state_ptr_ijk + ivar * state.nstride, 0,
                                     static_cast<int>(state.kstride)));
            });
        return d1;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    diff2(const amrex::Real *in_ptr, const int idx, const int stride) const
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
    mixed_diff2(const amrex::Real *in_ptr, const int idx, const int stride1,
                const int stride2) const
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

    /// Calculates all second derivatives for a single variable
    [[nodiscard]] AMREX_GPU_DEVICE
        AMREX_FORCE_INLINE amrex::Array2D<amrex::Real, 0, 3, 0, 3>
        diff2(int i, int j, int k,
              const amrex::Array4<amrex::Real const> &state,
              const int ivar) const
    {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
        amrex::Array2D<amrex::Real, 0, 3, 0, 3> d2;
        const auto *state_ptr_ijk = state.ptr(i, j, k);
        amrex::GpuArray<int, AMREX_SPACEDIM> strides{
            1, static_cast<int>(state.jstride),
            static_cast<int>(state.kstride)};
        const auto *pvar = state_ptr_ijk + ivar * state.nstride;
        FOR (dir1) // First calculate the repeated derivatives
        {
            d2(dir1, dir1) = diff2(pvar, 0, strides[dir1]);
            for (int dir2 = 0; dir2 < dir1; ++dir2)
            {
                auto tmp = mixed_diff2(pvar, 0, strides[dir1], strides[dir2]);
                d2(dir1, dir2) = tmp;
                d2(dir2, dir1) = tmp;
            }
        }
        return d2;
    }

    template <template <typename> class vars_t>
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE auto
    diff2(int i, int j, int k,
          const amrex::Array4<amrex::Real const> &state) const
    {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
        vars_t<amrex::Array2D<amrex::Real, 0, 3, 0, 3>> d2;
        const auto *state_ptr_ijk = state.ptr(i, j, k);
        amrex::GpuArray<int, AMREX_SPACEDIM> strides{
            1, static_cast<int>(state.jstride),
            static_cast<int>(state.kstride)};
        d2.enum_mapping(
            [&](const int &ivar, amrex::Array2D<amrex::Real, 0, 3, 0, 3> &var)
            {
                const auto *pvar = state_ptr_ijk + ivar * state.nstride;
                FOR (dir1) // First calculate the repeated derivatives
                {
                    var(dir1, dir1) = diff2(pvar, 0, strides[dir1]);
                    for (int dir2 = 0; dir2 < dir1; ++dir2)
                    {
                        auto tmp =
                            mixed_diff2(pvar, 0, strides[dir1], strides[dir2]);
                        var(dir1, dir2) = tmp;
                        var(dir2, dir1) = tmp;
                    }
                }
            });
        return d2;
    }

  protected: // Let's keep this protected ... we may want to change the
             // advection calculation
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    advection_term(const amrex::Real *in_ptr, const int idx,
                   const amrex::Real &vec_comp, const int stride,
                   const bool shift_positive) const
    {
        const amrex::Real in_left   = in_ptr[idx - stride];
        const amrex::Real in_centre = in_ptr[idx];
        const amrex::Real in_right  = in_ptr[idx + stride];

        amrex::Real weight_0 = -2.50000000000000000000e-1_rt;
        amrex::Real weight_1 = -8.33333333333333333333e-1_rt;
        amrex::Real weight_2 = +1.50000000000000000000e+0_rt;
        amrex::Real weight_3 = -5.00000000000000000000e-1_rt;
        amrex::Real weight_4 = +8.33333333333333333333e-2_rt;

        amrex::Real upwind =
            vec_comp *
            (weight_0 * in_left + weight_1 * in_centre + weight_2 * in_right +
             weight_3 * in_ptr[idx + 2 * stride] +
             weight_4 * in_ptr[idx + 3 * stride]) *
            m_one_over_dx;

        amrex::Real downwind =
            vec_comp *
            (-weight_4 * in_ptr[idx - 3 * stride] -
             weight_3 * in_ptr[idx - 2 * stride] - weight_2 * in_left -
             weight_1 * in_centre - weight_0 * in_right) *
            m_one_over_dx;

        return (shift_positive) ? upwind : downwind;
    }

  public:

    /// Calculates all second derivatives and returns as variable type specified
    /// by the template parameter
    template <template <typename> class vars_t>
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE auto
    advection(int i, int j, int k,
              const amrex::Array4<amrex::Real const> &state,
              const amrex::Array1D<amrex::Real, 0, 3> &vector) const
    {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
        vars_t<amrex::Real> advec;
        const auto *state_ptr_ijk = state.ptr(i, j, k);
        amrex::GpuArray<int, AMREX_SPACEDIM> strides{
            1, static_cast<int>(state.jstride),
            static_cast<int>(state.kstride)};
        advec.enum_mapping(
            [&](const int &ivar, amrex::Real &var)
            {
                var              = 0.;
                const auto *pvar = state_ptr_ijk + ivar * state.nstride;
                FOR (dir)
                {
                    const bool shift_positive = (vector(dir) > 0.0);
                    var += advection_term(pvar, 0, vector(dir), strides[dir],
                                          shift_positive);
                }
            });
        return advec;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    dissipation_term(const double *in_ptr, const int idx,
                     const int stride) const
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

    template <template <typename> class vars_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    add_dissipation(int i, int j, int k, vars_t<amrex::Real> &vars,
                    const amrex::Array4<amrex::Real const> &state,
                    const double factor) const
    {
        const auto *state_ptr_ijk = state.ptr(i, j, k);
        amrex::GpuArray<int, AMREX_SPACEDIM> strides{
            1, static_cast<int>(state.jstride),
            static_cast<int>(state.kstride)};
        vars.enum_mapping(
            [&](const int &ivar, amrex::Real &var)
            {
                FOR (dir)
                {
                    const auto stride  = strides[dir];
                    var               += factor *
                           dissipation_term(
                               state_ptr_ijk + ivar * state.nstride, 0, stride);
                }
            });
    }
};

#endif /* FOURTHORDERDERIVATIVES_HPP_ */
