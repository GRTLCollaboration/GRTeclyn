#ifndef LAGRANGEINTERPOLATION_HPP_
#define LAGRANGEINTERPOLATION_HPP_

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

template <int N> class LagrangeInterpolator
{

  private:
    int i0, j0,
        k0; // indices of the lower left corner of the stencil in the grid

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
    build_stencil(amrex::Real grid_pos, int &base_idx, amrex::Real *weights)
    {
        int stencil[N];
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
        amrex::Real rel_pos[N];
        for (int i = 0; i < N; ++i)
        {
            rel_pos[i] = static_cast<amrex::Real>(stencil[i]) - grid_pos;
        }

        // Compute Lagrange weights
        amrex::poly_interp_coeff(0.0, rel_pos, N, weights);

        // Write in the base (starting) index
        base_idx = stencil[0];
    }

  public:
    amrex::Real wx[N], wy[N],
        wz[N]; // where we store the weights for each dimension

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    LagrangeInterpolator() {};

    template <typename P>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
    compute_weights(const P &p,
                    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const &plo,
                    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const &dxi,
                    const amrex::IntVect &is_nodal)
    {

        // Compute the grid index of the position
        AMREX_D_TERM(
            amrex::Real lx =
                (amrex::Real(p.pos(0)) - plo[0]) * dxi[0] -
                static_cast<amrex::Real>(!is_nodal[0]) * amrex::Real(0.5);
            , amrex::Real ly =
                  (amrex::Real(p.pos(1)) - plo[1]) * dxi[1] -
                  static_cast<amrex::Real>(!is_nodal[1]) * amrex::Real(0.5);
            , amrex::Real lz =
                  (amrex::Real(p.pos(2)) - plo[2]) * dxi[2] -
                  static_cast<amrex::Real>(!is_nodal[2]) * amrex::Real(0.5););

        build_stencil(lx, i0, wx);

        // std::cout << "wx: " << wx[0] << ", " << wx[1] << ", " << wx[2] << ",
        // " << wx[3] << ", " << wx[4] << std::endl;

#if AMREX_SPACEDIM >= 2
        build_stencil(ly, j0, wy);
        ;
        // std::cout << "wy: " << wy[0] << ", " << wy[1] << ", " << wy[2] << ",
        // " << wy[3] << ", " << wy[4] << std::endl;
#endif
#if AMREX_SPACEDIM == 3
        build_stencil(lz, k0, wz);
        // std::cout << "wz: " << wz[0] << ", " << wz[1] << ", " << wz[2] << ",
        // " << wz[3] << ", " << wz[4] << std::endl;
#endif
    }

    // Function to perform the interpolation
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
    interpolate(const amrex::Array4<amrex::Real const> *data_arr,
                amrex::ParticleReal *val, int start_comp, int ncomp) const
    {
        int ctr          = 0;
        auto const &data = data_arr[0];

        for (int comp = start_comp; comp < start_comp + ncomp; ++comp)
        {
            val[ctr] = amrex::ParticleReal(0.0);
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
                        val[ctr] += data(amrex::IntVect(AMREX_D_DECL(
                                             i0 + ii, j0 + jj, k0 + kk)),
                                         comp) *
                                    AMREX_D_TERM(wx[ii], *wy[jj], *wz[kk]);
                    }
#if AMREX_SPACEDIM >= 2
                }
#endif
#if AMREX_SPACEDIM == 3
            }
#endif
            ++ctr;
        } // end of for comp loop

        // std::cout << "z at i = " << i0+0 << " j = " << j0+0 << " k = " <<
        // k0+0 << " : " << data(amrex::IntVect(AMREX_D_DECL(i0 + 0, j0 + 0, k0
        // + 0)), start_comp) << std::endl; std::cout << "z at i = " << i0+0 <<
        // " j = " << j0+0 << " k = " << k0+1 << " : " <<
        // data(amrex::IntVect(AMREX_D_DECL(i0 + 0, j0 + 0, k0 + 1)),
        // start_comp) << std::endl; std::cout << "z at i = " << i0+0 << " j = "
        // << j0+0 << " k = " << k0+2 << " : " <<
        // data(amrex::IntVect(AMREX_D_DECL(i0 + 0, j0 + 0, k0 + 2)),
        // start_comp) << std::endl; std::cout << "z at i = " << i0+0 << " j = "
        // << j0+0 << " k = " << k0+3 << " : " <<
        // data(amrex::IntVect(AMREX_D_DECL(i0 + 0, j0 + 0, k0 + 3)),
        // start_comp) << std::endl; std::cout << "z at i = " << i0+0 << " j = "
        // << j0+0 << " k = " << k0+4 << " : " <<
        // data(amrex::IntVect(AMREX_D_DECL(i0 + 0, j0 + 0, k0 + 4)),
        // start_comp) << std::endl;
    }
};

#endif /* LAGRANGEINTERPOLATION_HPP_ */
