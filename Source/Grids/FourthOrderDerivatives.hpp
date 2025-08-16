/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef FOURTHORDERDERIVATIVES_HPP_
#define FOURTHORDERDERIVATIVES_HPP_

#include "CCZ4Vars.hpp"
#include "Cell.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include <array>

class FourthOrderDerivatives
{
  private:
    amrex::Real m_dx;
    amrex::Real m_one_over_dx;
    amrex::Real m_one_over_dx2;

    // A function to return the right index for the tensor
    [[nodiscard]] AMREX_FORCE_INLINE const int var_idx(int ivar, int i,
                                                       int j) const
    {
        return ivar + i + j + ((i * j != 0) ? 1 : 0);
    }

  public:
    AMREX_GPU_HOST_DEVICE FourthOrderDerivatives(double dx)
        : m_dx(dx), m_one_over_dx(1 / dx), m_one_over_dx2(1 / (dx * dx))
    {
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    diff1(const amrex::Real *in_ptr, const int stride, const int idx = 0) const
    {
        amrex::Real weight_far  = 8.33333333333333333333e-2;
        amrex::Real weight_near = 6.66666666666666666667e-1;

        return (weight_far * in_ptr[idx - 2 * stride] -
                weight_near * in_ptr[idx - stride] +
                weight_near * in_ptr[idx + stride] -
                weight_far * in_ptr[idx + 2 * stride]) *
               m_one_over_dx;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor<1, amrex::Real>
    diff1(int ix, int iy, int iz, const amrex::Array4<const amrex::Real> &state,
          const int ivar) const
    {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
        Tensor<1, amrex::Real> d1;
        const auto *state_ptr_xyz = state.ptr(ix, iy, iz);
        amrex::GpuArray<int, AMREX_SPACEDIM> strides{
            1, static_cast<int>(state.jstride),
            static_cast<int>(state.kstride)};
        const auto *var_ptr = state_ptr_xyz + ivar * state.nstride;
        FOR (idir)
        {
            d1[idir] = diff1(var_ptr, strides[idir]);
        }
        return d1;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor<1, amrex::Real>
    diff1_scalar(int ix, int iy, int iz,
                 const amrex::Array4<const amrex::Real> &state,
                 const int ivar) const
    {
        return diff1(ix, iy, iz, state, ivar);
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor<2, amrex::Real>
    diff1_vector(int ix, int iy, int iz,
                 const amrex::Array4<const amrex::Real> &state,
                 const int ivar_0) const
    {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
        Tensor<2, amrex::Real> d1_vector;
        FOR (icomp)
        {
            const int ivar                = ivar_0 + icomp;
            Tensor<1, amrex::Real> d1_var = diff1(ix, iy, iz, state, ivar);
            FOR (idir)
            {
                d1_vector[icomp][idir] = d1_var[idir];
            }
        }
        return d1_vector;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor<3, amrex::Real>
    diff1_tensor(int ix, int iy, int iz,
                 const amrex::Array4<const amrex::Real> &state,
                 const int ivar_0) const
    {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
        Tensor<3, amrex::Real> d1_tensor;
        FOR (icomp, jcomp)
        {
            const int ivar                = var_idx(ivar_0, icomp, jcomp);
            Tensor<1, amrex::Real> d1_var = diff1(ix, iy, iz, state, ivar);
            FOR (idir)
            {
                d1_tensor[icomp][jcomp][idir] = d1_var[idir];
            }
        }
        return d1_tensor;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    diff2(const amrex::Real *in_ptr, const int stride, const int idx = 0) const
    {
        amrex::Real weight_far   = 8.33333333333333333333e-2;
        amrex::Real weight_near  = 1.33333333333333333333e+0;
        amrex::Real weight_local = 2.50000000000000000000e+0;

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
        amrex::Real weight_far_far   = 6.94444444444444444444e-3;
        amrex::Real weight_near_far  = 5.55555555555555555556e-2;
        amrex::Real weight_near_near = 4.44444444444444444444e-1;

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

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor<2, amrex::Real>
    diff2(int ix, int iy, int iz, const amrex::Array4<amrex::Real const> &state,
          const int ivar) const
    {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
        Tensor<2, amrex::Real> d2;
        const auto *state_ptr_xyz = state.ptr(ix, iy, iz);
        amrex::GpuArray<int, AMREX_SPACEDIM> strides{
            1, static_cast<int>(state.jstride),
            static_cast<int>(state.kstride)};
        const auto *var_ptr = state_ptr_xyz + ivar * state.nstride;
        FOR (idir1)
        {
            d2[idir1][idir1] = diff2(var_ptr, strides[idir1]);
            for (int idir2 = 0; idir2 < idir1; ++idir2)
            {
                auto d2_tmp =
                    mixed_diff2(var_ptr, strides[idir1], strides[idir2]);
                d2[idir1][idir2] = d2_tmp;
                d2[idir2][idir1] = d2_tmp;
            }
        }
        return d2;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor<2, amrex::Real>
    diff2_scalar(int ix, int iy, int iz,
                 const amrex::Array4<const amrex::Real> &state,
                 const int ivar) const
    {
        return diff2(ix, iy, iz, state, ivar);
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor<3, amrex::Real>
    diff2_vector(int ix, int iy, int iz,
                 const amrex::Array4<const amrex::Real> &state,
                 const int ivar_0) const
    {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
        Tensor<3, amrex::Real> d2_vector;
        FOR (icomp)
        {
            const int ivar                = ivar_0 + icomp;
            Tensor<2, amrex::Real> d2_var = diff2(ix, iy, iz, state, ivar);
            FOR (idir, jdir)
            {
                d2_vector[icomp][idir][jdir] = d2_var[idir][jdir];
            }
        }
        return d2_vector;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor<4, amrex::Real>
    diff2_tensor(int ix, int iy, int iz,
                 const amrex::Array4<const amrex::Real> &state,
                 const int ivar_0) const
    {
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
        Tensor<4, amrex::Real> d2_tensor;
        FOR (icomp, jcomp)
        {
            const int ivar                = var_idx(ivar_0, icomp, jcomp);
            Tensor<2, amrex::Real> d1_var = diff2(ix, iy, iz, state, ivar);

            FOR (idir, jdir)
            {
                d2_tensor[icomp][jcomp][idir][jdir] = d1_var[idir][jdir];
            }
        }
        return d2_tensor;
    }

  protected:
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    advection_term(const amrex::Real *in_ptr, const amrex::Real &shift_comp,
                   const int stride, const bool shift_positive,
                   const int idx = 0) const
    {
        amrex::Real weight_0 = -2.50000000000000000000e-1;
        amrex::Real weight_1 = -8.33333333333333333333e-1;
        amrex::Real weight_2 = +1.50000000000000000000e+0;
        amrex::Real weight_3 = -5.00000000000000000000e-1;
        amrex::Real weight_4 = +8.33333333333333333333e-2;

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
              const Tensor<1, amrex::Real> &shift_vector, const int ivar) const
    {
        amrex::Real advec         = 0.0;
        const auto *state_ptr_xyz = state.ptr(ix, iy, iz);
        amrex::GpuArray<int, AMREX_SPACEDIM> strides{
            1, static_cast<int>(state.jstride),
            static_cast<int>(state.kstride)};

        const auto *var_ptr = state_ptr_xyz + ivar * state.nstride;
        FOR (idir)
        {
            const bool shift_positive = (shift_vector[idir] > 0.0);
            advec += advection_term(var_ptr, shift_vector[idir], strides[idir],
                                    shift_positive);
        }
        return advec;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real advec_scalar(
        int ix, int iy, int iz, const amrex::Array4<amrex::Real const> &state,
        const Tensor<1, amrex::Real> &shift_vector, const int ivar) const
    {
        return advection(ix, iy, iz, state, shift_vector, ivar);
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor<1, amrex::Real>
    advec_vector(int ix, int iy, int iz,
                 const amrex::Array4<amrex::Real const> &state,
                 const Tensor<1, amrex::Real> &shift_vector,
                 const int ivar0) const
    {
        Tensor<1, amrex::Real> advec_vector;
        FOR (icomp)
        {
            int ivar = ivar0 + icomp;
            advec_vector[icomp] =
                advection(ix, iy, iz, state, shift_vector, ivar);
        }
        return advec_vector;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor<2, amrex::Real>
    advec_tensor(int ix, int iy, int iz,
                 const amrex::Array4<amrex::Real const> &state,
                 const Tensor<1, amrex::Real> &shift_vector,
                 const int ivar0) const
    {
        Tensor<2, amrex::Real> advec_tensor;
        FOR (icomp, jcomp)
        {
            int ivar = var_idx(ivar0, icomp, jcomp);
            advec_tensor[icomp][jcomp] =
                advection(ix, iy, iz, state, shift_vector, ivar);
        }
        return advec_tensor;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    dissipation_term(const double *in_ptr, const int stride,
                     const int idx = 0) const
    {
        amrex::Real weight_vfar  = 1.56250e-2;
        amrex::Real weight_far   = 9.37500e-2;
        amrex::Real weight_near  = 2.34375e-1;
        amrex::Real weight_local = 3.12500e-1;

        return (weight_vfar * in_ptr[idx - 3 * stride] -
                weight_far * in_ptr[idx - 2 * stride] +
                weight_near * in_ptr[idx - stride] -
                weight_local * in_ptr[idx] +
                weight_near * in_ptr[idx + stride] -
                weight_far * in_ptr[idx + 2 * stride] +
                weight_vfar * in_ptr[idx + 3 * stride]) *
               m_one_over_dx;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    calculate_dissipation(int ix, int iy, int iz,
                          const amrex::Array4<amrex::Real const> &state,
                          const double sigma_coeff, const int ivar) const
    {
        amrex::Real diss          = 0.0;
        const auto *state_ptr_xyz = state.ptr(ix, iy, iz);
        amrex::GpuArray<int, AMREX_SPACEDIM> strides{
            1, static_cast<int>(state.jstride),
            static_cast<int>(state.kstride)};
        FOR (idir)
        {
            const auto stride = strides[idir];
            diss +=
                sigma_coeff *
                dissipation_term(state_ptr_xyz + ivar * state.nstride, stride);
        }
        return diss;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    add_dissipation(int ix, int iy, int iz, CCZ4Vars &vars,
                    const amrex::Array4<amrex::Real const> &state,
                    const double sigma_coeff,
                    int num_vars = NUM_CCZ4_VARS) const
    {
        for (int ivar = 0; ivar < num_vars; ivar++)
        {
            amrex::Real diss =
                calculate_dissipation(ix, iy, iz, state, sigma_coeff, ivar);
            amrex::Real var_plus_diss = vars.get_var(ivar) + diss;
            vars.store_var(ivar, var_plus_diss);
        }
    }
};

#endif /* FOURTHORDERDERIVATIVES2_HPP_ */
