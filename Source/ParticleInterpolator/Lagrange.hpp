/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef LAGRANGE_HPP_
#define LAGRANGE_HPP_

#include "InterpolationQueryParticle.hpp"
#include <AMReX_Array4.H>
#include <AMReX_Gpu.H>
#include <AMReX_IntVect.H>
#include <AMReX_REAL.H>
#include <cmath>
#include "AMReX_LOUtil_K.H"

#include <limits>

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
    build_stencil(amrex::Real grid_pos, int &base_idx, amrex::Real *weights,
                  bool lo_reflective, bool hi_reflective, int ncell,
                  int deriv_order)
    {
        // Build interpolation stencil for a particular derivative:
        // 0 = no derivative
        // 1 = 1st derivative
        // 2 = 2nd derivative
        std::array<int, N> stencil;
        constexpr int center_offset =
            N / 2; // offset based on the number of points we are using
        int center{};
        constexpr amrex::Real eps =
            10 * std::numeric_limits<
                     amrex::Real>::epsilon(); // choose a small number around
                                              // machine round-off precision
        // when the position is very close to the axis of symmetry and
        // reflective boundary conditions are used, we can run into some
        // problems of not having enough ghosts. For example, for low symmetric
        // boundary, the center of the stencil can end up being rounded up to
        // -1. If center = -1, then the stencil is e.g. (-3, -2, -1, 0, 1) for
        // 4th order interpolation, which is problematic here as -3 is out of
        // bounds (note that we fill [4/2]=2 ghost cells). To avoid this, we
        // will default to center = 0 in this situation.

        const auto lo_face = amrex::Real(
            -0.5); // for low symmetric boundary, the face is at -0.5
        const auto hi_face =
            amrex::Real(ncell) -
            amrex::Real(
                0.5); // for high symmetric boundary, the face is at ncell-0.5

        if (lo_reflective && amrex::Math::abs(grid_pos - lo_face) < eps)
        {
            center = 0;
        }
        else if (hi_reflective && amrex::Math::abs(grid_pos - hi_face) < eps)
        {
            center = ncell - 1;
        }
        else
        {
            center = static_cast<int>(
                amrex::Math::round(grid_pos)); // round the grid position to get
                                               // the nearest cell center
        }

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

        // Compute Lagrange weights for the requested derivative order
        switch (deriv_order)
        {
        case 1:
            poly_interp_coeff_d1(0.0, rel_pos, N, weights);
            break;
        case 2:
            poly_interp_coeff_d2(0.0, rel_pos, N, weights);
            break;
        default:
            amrex::poly_interp_coeff(0.0, rel_pos, N, weights);
            break;
        }

        // Write in the base (starting) index
        base_idx = stencil[0];
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
    poly_interp_coeff_d1(amrex::Real xInt, amrex::Real const *AMREX_RESTRICT x,
                         int order,
                         amrex::Real *AMREX_RESTRICT weights) noexcept
    {
        for (int j = 0; j < order; ++j)
        {
            auto den = amrex::Real(1.0);
            for (int i = 0; i < order; ++i)
            {
                if (i != j)
                {
                    den *= (x[j] - x[i]);
                }
            }

            weights[j] = amrex::Real(0.0);

            for (int k = 0; k < order; ++k)
            {
                if (k == j)
                {
                    continue;
                }
                auto num = amrex::Real(1.0);
                for (int i = 0; i < order; ++i)
                {
                    if (i != j && i != k)
                    {
                        num *= (xInt - x[i]);
                    }
                }

                weights[j] += num / den;
            }
        }
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
    poly_interp_coeff_d2(amrex::Real xInt, amrex::Real const *AMREX_RESTRICT x,
                         int order,
                         amrex::Real *AMREX_RESTRICT weights) noexcept
    {
        for (int j = 0; j < order; ++j)
        {
            weights[j] = amrex::Real(0.0);

            auto den = amrex::Real(1.0);
            for (int i = 0; i < order; ++i)
            {
                if (i != j)
                {
                    den *= (x[j] - x[i]);
                }
            }

            for (int l = 0; l < order; ++l)
            {
                if (l == j)
                {
                    continue;
                }

                for (int k = 0; k < order; ++k)
                {
                    if (k == j)
                    {
                        continue;
                    }
                    if (k == l)
                    {
                        continue;
                    }

                    auto num = amrex::Real(1.0);

                    for (int i = 0; i < order; ++i)
                    {
                        if (i != j && i != k && i != l)
                        {
                            num *= (xInt - x[i]);
                        }
                    }

                    weights[j] += num / den;
                }
            }
        }
    }

  public:
    // where we store the weights for each dimension
    // NOLINTBEGIN(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
    amrex::Real weights_local[AMREX_SPACEDIM][N]{};
    amrex::Real weights_d1[AMREX_SPACEDIM][N]{};
    amrex::Real weights_d2[AMREX_SPACEDIM][N]{};
    // NOLINTEND(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE Lagrange() = default;

    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    template <typename P>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    compute_weights(const P &par,
                    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const &plo,
                    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const &dxi,
                    const amrex::IntVect &is_nodal,
                    amrex::GpuArray<bool, AMREX_SPACEDIM> const lo_reflective,
                    amrex::GpuArray<bool, AMREX_SPACEDIM> const hi_reflective,
                    amrex::GpuArray<int, AMREX_SPACEDIM> const domain_ncell)
    // NOLINTEND(bugprone-easily-swappable-parameters)
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

        build_stencil(xpos, i0, weights_local[0], lo_reflective[1],
                      hi_reflective[1], domain_ncell[1], 0);
        build_stencil(xpos, i0, weights_d1[0], lo_reflective[1],
                      hi_reflective[1], domain_ncell[1], 1);
        build_stencil(xpos, i0, weights_d2[0], lo_reflective[1],
                      hi_reflective[1], domain_ncell[1], 2);
#if AMREX_SPACEDIM >= 2
        build_stencil(ypos, j0, weights_local[1], lo_reflective[1],
                      hi_reflective[1], domain_ncell[1], 0);
        build_stencil(ypos, j0, weights_d1[1], lo_reflective[1],
                      hi_reflective[1], domain_ncell[1], 1);
        build_stencil(ypos, j0, weights_d2[1], lo_reflective[1],
                      hi_reflective[1], domain_ncell[1], 2);

#endif
#if AMREX_SPACEDIM == 3
        build_stencil(zpos, k0, weights_local[2], lo_reflective[1],
                      hi_reflective[1], domain_ncell[1], 0);
        build_stencil(zpos, k0, weights_d1[2], lo_reflective[1],
                      hi_reflective[1], domain_ncell[1], 1);
        build_stencil(zpos, k0, weights_d2[2], lo_reflective[1],
                      hi_reflective[1], domain_ncell[1], 2);
#endif
    }

    // Function to perform the interpolation
    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    interpolate(const amrex::Array4<amrex::Real const> *data_arr,
                amrex::ParticleReal *val, const Derivative *derivs,
                InterpolationQueryParticle::out_t *const *comps_arr,
                const int *comp_counts, int ncomp, amrex::Real const dx) const
    // NOLINTEND(bugprone-easily-swappable-parameters)
    {

        int counter      = 0;
        auto const &data = data_arr[0];
        // NOLINTBEGIN(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
        amrex::Real weights[AMREX_SPACEDIM][N]{};
        // NOLINTEND(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)

        for (int i = 0; i < ncomp; ++i)
        {
            const Derivative &deriv = derivs[i];

            const InterpolationQueryParticle::out_t *comps = comps_arr[i];

            // Collect dimensions where we want at least a first derivative
            for (int dim = 0; dim < AMREX_SPACEDIM; ++dim)
            {
                if (deriv[dim] == 1)
                {
                    for (int j = 0; j < N; j++)
                    {
                        weights[dim][j] = weights_d1[dim][j] / dx;
                    }
                }
                else if (deriv[dim] == 2)
                {
                    for (int j = 0; j < N; j++)
                    {
                        weights[dim][j] = weights_d2[dim][j] / pow(dx, 2);
                    }
                }
                else
                {
                    for (int j = 0; j < N; j++)
                    {
                        weights[dim][j] = weights_local[dim][j];
                    }
                }
            }

            for (int j = 0; j < comp_counts[i]; ++j)
            {
                const int comp = comps[j].comp;

                val[counter] = amrex::ParticleReal(0.0);

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
                            val[counter] +=
                                data(amrex::IntVect(AMREX_D_DECL(
                                         i0 + ii, j0 + jj, k0 + kk)),
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
                ++counter;
            }
        }
    }
};
#endif /* LAGRANGE_HPP_ */
