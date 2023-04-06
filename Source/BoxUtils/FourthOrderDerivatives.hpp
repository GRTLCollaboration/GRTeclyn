/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FOURTHORDERDERIVATIVES_HPP_
#define FOURTHORDERDERIVATIVES_HPP_

#include "Cell.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp"
#include "simd.hpp"
#include <array>

class FourthOrderDerivatives
{
  private:
    double m_dx;
    double m_one_over_dx;
    double m_one_over_dx2;

  public:
    FourthOrderDerivatives(double dx)
        : m_dx(dx), m_one_over_dx(1 / dx), m_one_over_dx2(1 / (dx * dx))
    {
    }

    template <class data_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE data_t diff1(const double *in_ptr,
                                                     const int idx,
                                                     const int stride) const
    {
        const auto *in = SIMDIFY<data_t>(in_ptr);

        data_t weight_far  = 8.33333333333333333333e-2;
        data_t weight_near = 6.66666666666666666667e-1;

        // NOTE: if you have been sent here by the debugger because of
        // EXC_BAD_ACCESS  or something similar you might be trying to take
        // derivatives without ghost points.
        return (weight_far * in[idx - 2 * stride] -
                weight_near * in[idx - stride] +
                weight_near * in[idx + stride] -
                weight_far * in[idx + 2 * stride]) *
               m_one_over_dx;
    }

    /// Calculates all first derivatives and returns as variable type specified
    /// by the template parameter
    template <template <typename> class vars_t, class data_t>
    AMREX_GPU_DEVICE [[nodiscard]] AMREX_FORCE_INLINE auto
    diff1(int i, int j, int k, const amrex::Array4<data_t const> &state) const
    {
        vars_t<Tensor<1, data_t>> d1;
        const auto *state_ptr_ijk = state.ptr(i, j, k);
        d1.enum_mapping(
            [&](const int &ivar, Tensor<1, data_t> &var)
            {
                AMREX_D_TERM(
                    var[0] = diff1<data_t>(state_ptr_ijk + ivar * state.nstride,
                                           0, 1);
                    ,
                    var[1] = diff1<data_t>(state_ptr_ijk + ivar * state.nstride,
                                           0, static_cast<int>(state.jstride));
                    ,
                    var[2] = diff1<data_t>(state_ptr_ijk + ivar * state.nstride,
                                           0, static_cast<int>(state.kstride)));
            });
        return d1;
    }

    template <class data_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE data_t diff2(const double *in_ptr,
                                                     const int idx,
                                                     const int stride) const
    {
        const auto *in = SIMDIFY<data_t>(in_ptr);

        data_t weight_far   = 8.33333333333333333333e-2;
        data_t weight_near  = 1.33333333333333333333e+0;
        data_t weight_local = 2.50000000000000000000e+0;

        return (-weight_far * in[idx - 2 * stride] +
                weight_near * in[idx - stride] - weight_local * in[idx] +
                weight_near * in[idx + stride] -
                weight_far * in[idx + 2 * stride]) *
               m_one_over_dx2;
    }

    template <class data_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE data_t
    mixed_diff2(const double *in_ptr, const int idx, const int stride1,
                const int stride2) const
    {
        const auto *in = SIMDIFY<data_t>(in_ptr);

        data_t weight_far_far   = 6.94444444444444444444e-3;
        data_t weight_near_far  = 5.55555555555555555556e-2;
        data_t weight_near_near = 4.44444444444444444444e-1;

        return (weight_far_far * in[idx - 2 * stride1 - 2 * stride2] -
                weight_near_far * in[idx - 2 * stride1 - stride2] +
                weight_near_far * in[idx - 2 * stride1 + stride2] -
                weight_far_far * in[idx - 2 * stride1 + 2 * stride2]

                - weight_near_far * in[idx - stride1 - 2 * stride2] +
                weight_near_near * in[idx - stride1 - stride2] -
                weight_near_near * in[idx - stride1 + stride2] +
                weight_near_far * in[idx - stride1 + 2 * stride2]

                + weight_near_far * in[idx + stride1 - 2 * stride2] -
                weight_near_near * in[idx + stride1 - stride2] +
                weight_near_near * in[idx + stride1 + stride2] -
                weight_near_far * in[idx + stride1 + 2 * stride2]

                - weight_far_far * in[idx + 2 * stride1 - 2 * stride2] +
                weight_near_far * in[idx + 2 * stride1 - stride2] -
                weight_near_far * in[idx + 2 * stride1 + stride2] +
                weight_far_far * in[idx + 2 * stride1 + 2 * stride2]) *
               m_one_over_dx2;
    }

    /// Calculates all second derivatives and returns as variable type specified
    /// by the template parameter
    template <template <typename> class vars_t, class data_t>
    AMREX_GPU_DEVICE [[nodiscard]] AMREX_FORCE_INLINE auto
    diff2(int i, int j, int k, amrex::Array4<data_t const> const &state) const
    {
        vars_t<Tensor<2, data_t>> d2{};
        const auto *state_ptr_ijk = state.ptr(i, j, k);
        amrex::GpuArray<int, AMREX_SPACEDIM> strides{
            1, static_cast<int>(state.jstride),
            static_cast<int>(state.kstride)};
        d2.enum_mapping(
            [&](const int &ivar, Tensor<2, data_t> &var)
            {
                const auto *pvar = state_ptr_ijk + ivar * state.nstride;
                FOR (dir1) // First calculate the repeated derivatives
                {
                    var[dir1][dir1] = diff2<data_t>(pvar, 0, strides[dir1]);
                    for (int dir2 = 0; dir2 < dir1; ++dir2)
                    {
                        auto tmp = mixed_diff2<data_t>(pvar, 0, strides[dir1],
                                                       strides[dir2]);
                        var[dir1][dir2] = tmp;
                        var[dir2][dir1] = tmp;
                    }
                }
            });
        return d2;
    }

  protected: // Let's keep this protected ... we may want to change the
             // advection calculation
    template <class data_t, class mask_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE data_t
    advection_term(const double *in_ptr, const int idx, const data_t &vec_comp,
                   const int stride, const mask_t shift_positive) const
    {
        const auto *const in   = SIMDIFY<data_t>(in_ptr);
        const data_t in_left   = in[idx - stride];
        const data_t in_centre = in[idx];
        const data_t in_right  = in[idx + stride];

        data_t weight_0 = -2.50000000000000000000e-1;
        data_t weight_1 = -8.33333333333333333333e-1;
        data_t weight_2 = +1.50000000000000000000e+0;
        data_t weight_3 = -5.00000000000000000000e-1;
        data_t weight_4 = +8.33333333333333333333e-2;

        data_t upwind;
        upwind = vec_comp *
                 (weight_0 * in_left + weight_1 * in_centre +
                  weight_2 * in_right + weight_3 * in[idx + 2 * stride] +
                  weight_4 * in[idx + 3 * stride]) *
                 m_one_over_dx;

        data_t downwind;
        downwind = vec_comp *
                   (-weight_4 * in[idx - 3 * stride] -
                    weight_3 * in[idx - 2 * stride] - weight_2 * in_left -
                    weight_1 * in_centre - weight_0 * in_right) *
                   m_one_over_dx;

        return simd_conditional(shift_positive, upwind, downwind);
    }

  public:

    /// Calculates all second derivatives and returns as variable type specified
    /// by the template parameter
    template <template <typename> class vars_t, class data_t>
    AMREX_GPU_DEVICE [[nodiscard]] AMREX_FORCE_INLINE auto
    advection(int i, int j, int k, amrex::Array4<data_t const> const &state,
              const Tensor<1, data_t> &vector) const
    {
        vars_t<data_t> advec;
        const auto *state_ptr_ijk = state.ptr(i, j, k);
        amrex::GpuArray<int, AMREX_SPACEDIM> strides{
            1, static_cast<int>(state.jstride),
            static_cast<int>(state.kstride)};
        advec.enum_mapping(
            [&](const int &ivar, data_t &var)
            {
                var              = 0.;
                const auto *pvar = state_ptr_ijk + ivar * state.nstride;
                FOR (dir)
                {
                    const auto shift_positive =
                        simd_compare_gt(vector[dir], 0.0);
                    var += advection_term(pvar, 0, vector[dir], strides[dir],
                                          shift_positive);
                }
            });
        return advec;
    }

    template <class data_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE data_t dissipation_term(
        const double *in_ptr, const int idx, const int stride) const
    {
        const auto *const in = SIMDIFY<data_t>(in_ptr);
        data_t weight_vfar   = 1.56250e-2;
        data_t weight_far    = 9.37500e-2;
        data_t weight_near   = 2.34375e-1;
        data_t weight_local  = 3.12500e-1;

        return (weight_vfar * in[idx - 3 * stride] -
                weight_far * in[idx - 2 * stride] +
                weight_near * in[idx - stride] - weight_local * in[idx] +
                weight_near * in[idx + stride] -
                weight_far * in[idx + 2 * stride] +
                weight_vfar * in[idx + 3 * stride]) *
               m_one_over_dx;
    }

    template <class data_t, template <typename> class vars_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    add_dissipation(int i, int j, int k, vars_t<data_t> &vars,
                    amrex::Array4<data_t const> const &state,
                    const double factor) const
    {
        const auto *state_ptr_ijk = state.ptr(i, j, k);
        amrex::GpuArray<int, AMREX_SPACEDIM> strides{
            1, static_cast<int>(state.jstride),
            static_cast<int>(state.kstride)};
        vars.enum_mapping(
            [&](const int &ivar, data_t &var)
            {
                FOR (dir)
                {
                    const auto stride = strides[dir];
                    var               += factor *
                           dissipation_term<data_t>(
                               state_ptr_ijk + ivar * state.nstride, 0, stride);
                }
            });
    }
};

#endif /* FOURTHORDERDERIVATIVES_HPP_ */
