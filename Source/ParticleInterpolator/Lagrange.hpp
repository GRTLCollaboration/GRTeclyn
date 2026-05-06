/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef LAGRANGE_HPP_
#define LAGRANGE_HPP_

#include "FourthOrderDerivatives.hpp"
#include "InterpolationQueryParticle.hpp"
#include <AMReX_Array4.H>
#include <AMReX_Gpu.H>
#include <AMReX_IntVect.H>
#include <AMReX_REAL.H>
#include <cmath>
#include "AMReX_LOUtil_K.H"

// Class for (N-1)th order interpolation of the mesh data onto the particle
// using Lagrange polynomials. Currently, it allows to interpolate only one
// field at a time. But the field itself can have multiple components. Assumes
// uniform grids

template <int N> class Lagrange
{
  private:
    // indices of the lower left corner of the stencil in the grid
    int i0{};
    int j0{};
    int k0{};

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    build_stencil(amrex::Real grid_pos, int &base_idx, amrex::Real *weights)
    {
        std::array<int, N> stencil;
        constexpr int center_offset =
            N / 2; // offset based on the number of points we are using
        int center = static_cast<int>(
            amrex::Math::round(grid_pos)); // round the grid position to get the
                                           // nearest cell center

        // Fill in the stencil around the position of interpolation
        for (int i = 0; i < N; ++i)
        {
            stencil[i] = center - center_offset + i;
        }

        // Compute the relative position from the interpolation point
        // NOLINTBEGIN(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
        amrex::Real rel_pos[N];
        // NOLINTEND(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
        for (int i = 0; i < N; ++i)
        {
            rel_pos[i] = static_cast<amrex::Real>(stencil[i]) - grid_pos;
        }

        // Compute Lagrange weights
        amrex::poly_interp_coeff(0.0, rel_pos, N, weights);

        // Write in the base (starting) index
        base_idx = stencil[0];
    }

    // Different cases for interpolation of derivatives
    amrex::ParticleReal
    interp_deriv_local(int comp,
                       const amrex::Array4<amrex::Real const> data) const
    {
        amrex::ParticleReal val = amrex::ParticleReal(0.0);

#if AMREX_SPACEDIM == 3
        for (int kk = 0; kk < N; ++kk)
        {
#endif
#if AMREX_SPACEDIM >= 2
            for (int jj = 0; jj < N; ++jj)
            {
#endif
                for (int ii = 0; ii < N; ++ii)
                {
                    val += data(amrex::IntVect(
                                    AMREX_D_DECL(i0 + ii, j0 + jj, k0 + kk)),
                                comp) *
                           AMREX_D_TERM(weights[0][ii], *weights[1][jj],
                                        *weights[2][kk]);
                }
#if AMREX_SPACEDIM >= 2
            }
#endif
#if AMREX_SPACEDIM == 3
        }
#endif
        return val;
    }

    amrex::ParticleReal
    interp_deriv_1d(int comp, int dim, int stride,
                    FourthOrderDerivatives a_deriv,
                    const amrex::Array4<amrex::Real const> data) const
    {
        amrex::ParticleReal val = amrex::ParticleReal(0.0);

        for (int ii = 0; ii < N; ++ii)
        {
            val += a_deriv.diff1(data.ptr(i0, j0, k0) +
                                     comp * data.stride.a[2] + ii * stride,
                                 stride) *
                   weights[dim][ii];
        }
        return val;
    }

    amrex::ParticleReal
    interp_deriv_2d(int comp, int dim, int stride,
                    FourthOrderDerivatives a_deriv,
                    const amrex::Array4<amrex::Real const> data) const
    {
        amrex::ParticleReal val = amrex::ParticleReal(0.0);

        for (int ii = 0; ii < N; ++ii)
        {
            val += a_deriv.diff2(data.ptr(i0, j0, k0) +
                                     comp * data.stride.a[2] + ii * stride,
                                 stride) *
                   weights[dim][ii];
        }
        return val;
    }

    amrex::ParticleReal
    interp_deriv_2d_mixed(int comp, int dim1, int dim2, int stride1,
                          int stride2, FourthOrderDerivatives a_deriv,
                          const amrex::Array4<amrex::Real const> data) const
    {
        amrex::ParticleReal val = amrex::ParticleReal(0.0);

        for (int ii = 0; ii < N; ++ii)
        {
            for (int jj = 0; jj < N; ++jj)
            {
                val += a_deriv.mixed_diff2(data.ptr(i0, j0, k0) +
                                               comp * data.stride.a[2] +
                                               ii * stride1 + jj * stride2,
                                           stride1, stride2) *
                       weights[dim1][ii] * weights[dim2][jj];
            }
        }
        return val;
    }

  public:
    // where we store the weights for each dimension
    // NOLINTBEGIN(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
    amrex::Real weights[AMREX_SPACEDIM][N]{};
    // NOLINTEND(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE Lagrange() = default;

    template <typename P>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    compute_weights(const P &par,
                    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const &plo,
                    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const &dxi,
                    const amrex::IntVect &is_nodal)
    {

        // Compute the grid index of the position
        AMREX_D_TERM(
            amrex::Real xpos =
                (amrex::Real(par.pos(0)) - plo[0]) * dxi[0] -
                static_cast<amrex::Real>(!is_nodal[0]) * amrex::Real(0.5);
            , amrex::Real ypos =
                  (amrex::Real(par.pos(1)) - plo[1]) * dxi[1] -
                  static_cast<amrex::Real>(!is_nodal[1]) * amrex::Real(0.5);
            , amrex::Real zpos =
                  (amrex::Real(par.pos(2)) - plo[2]) * dxi[2] -
                  static_cast<amrex::Real>(!is_nodal[2]) * amrex::Real(0.5););

        build_stencil(xpos, i0, weights[0]);
#if AMREX_SPACEDIM >= 2
        build_stencil(ypos, j0, weights[1]);
#endif
#if AMREX_SPACEDIM == 3
        build_stencil(zpos, k0, weights[2]);
#endif
    }

    // Function to perform the interpolation
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    interpolate(const amrex::Array4<amrex::Real const> *data_arr,
                amrex::ParticleReal *val,
                InterpolationQueryParticle::iterator deriv_it_begin,
                InterpolationQueryParticle::iterator deriv_it_end,
                amrex::Real const dx) const
    {
        int counter      = 0;
        auto const &data = data_arr[0];

        for (auto deriv_it = deriv_it_begin; deriv_it != deriv_it_end;
             ++deriv_it)
        {
            using comps_t =
                std::vector<typename InterpolationQueryParticle::out_t>;

            const Derivative &deriv = deriv_it->first;

            amrex::GpuArray<int, AMREX_SPACEDIM> strides{1, data.stride.a[0],
                                                         data.stride.a[1]};

            FourthOrderDerivatives a_deriv(dx);

            std::function<amrex::ParticleReal(int)> interp_deriv;

            comps_t &comps = deriv_it->second;

            // Select which interpolation function to call based on derivative
            if (deriv == Derivative::LOCAL)
            {
                interp_deriv = [=](int comp) -> amrex::ParticleReal
                { return interp_deriv_local(comp, data); };
            }
            else
            {
                std::vector<int> deriv_dims;

                // Collect dimensions where we want at least a first derivative
                for (int dim = 0; dim < AMREX_SPACEDIM; ++dim)
                {
                    if (deriv[dim] > 0)
                    {
                        deriv_dims.push_back(dim);
                    }
                }

                if (deriv_dims.size() == 1)
                {
                    // 2nd derivative in 1 dimension
                    if (deriv[deriv_dims[0]] == 2)
                    {
                        interp_deriv = [=](int comp) -> amrex::ParticleReal
                        {
                            return interp_deriv_2d(comp, deriv_dims[0],
                                                   strides[deriv_dims[0]],
                                                   a_deriv, data);
                        };
                    }
                    // 1st derivative
                    else
                    {
                        interp_deriv = [=](int comp) -> amrex::ParticleReal
                        {
                            return interp_deriv_1d(comp, deriv_dims[0],
                                                   strides[deriv_dims[0]],
                                                   a_deriv, data);
                        };
                    }
                }
                // Mixed 2nd derivative
                else if (deriv_dims.size() == 2)
                {
                    interp_deriv = [=](int comp) -> amrex::ParticleReal
                    {
                        return interp_deriv_2d_mixed(
                            comp, deriv_dims[0], deriv_dims[1],
                            strides[deriv_dims[0]], strides[deriv_dims[1]],
                            a_deriv, data);
                    };
                }
            }

            // Interpolate for comps with this derivative
            for (auto &entry : comps)
            {
                const int comp = entry.comp;

                val[counter] = interp_deriv(comp);
            }

            ++counter;
        }
    }
};
#endif /* LAGRANGE_HPP_ */
