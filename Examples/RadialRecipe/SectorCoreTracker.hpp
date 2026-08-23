/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SECTORCORETRACKER_HPP_
#define SECTORCORETRACKER_HPP_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Reduce.H>

#include <array>
#include <cmath>
#include <limits>

/// Halo-free position of one matter sector.
///
/// SINGLE RESPONSIBILITY.  This file answers exactly one question -- "where is
/// this sector's core?" -- and knows nothing about why anyone is asking.  It
/// does not shift grids, does not write files and does not read parameters; the
/// caller hands it the component indices it should weight by, so it is not tied
/// to any particular matter layout (DEPENDENCY INVERSION: the physics module
/// depends on this small interface, not the other way round).
///
/// The core definition deliberately matches the one the plotfile consumer uses
/// for `sector_dynamics.dat`
/// (`grteclyn-wrapper/.../extraction/sector_dynamics.py`):
///
///     w      = sqrt(sum of squares of the sector's four field components)
///     core   = cells with w >= core_fraction * max(w)
///     centre = centroid of the core weighted by w^4
///
/// so the in-situ number and the post-hoc number are the same quantity computed
/// the same way, and disagreement between them is a real signal rather than a
/// difference of definition.  Everything below works with s = w^2 to avoid a
/// square root per cell: `w >= f*w_max` is `s >= f^2*s_max`, and `w^4` is `s^2`.
struct SectorFieldSet
{
    //! Components whose squares sum to s = w^2.  Unused slots are set to -1.
    std::array<int, 4> comps{-1, -1, -1, -1};
};

/// Where a sector's core is, and how sure we are that there is one.
struct SectorCore
{
    double x{0.0};
    double y{0.0};
    double z{0.0};
    double peak{0.0};  //!< max w over the search region
    bool found{false}; //!< false => no matter above threshold; hold last value
};

/// Locate one sector's core on a single-level state.
///
/// `a_seed_*` and `a_ball_radius` restrict the search to a ball around the
/// previously known position, which keeps departed wake from pulling the answer
/// backwards.  Pass `a_ball_radius <= 0` to search the whole domain -- which is
/// what the first call should do, since there is no previous position yet.
inline SectorCore
track_sector_core(const amrex::MultiFab &a_state, const amrex::Geometry &a_geom,
                  const SectorFieldSet &a_fields, double a_seed_x,
                  double a_seed_y, double a_seed_z, double a_ball_radius,
                  double a_core_fraction = 0.25)
{
    const auto prob_lo = a_geom.ProbLoArray();
    const auto dx_arr  = a_geom.CellSizeArray();

    const amrex::Real ball_r2 =
        (a_ball_radius > 0.0)
            ? static_cast<amrex::Real>(a_ball_radius * a_ball_radius)
            : std::numeric_limits<amrex::Real>::max();
    const amrex::Real seed_x = static_cast<amrex::Real>(a_seed_x);
    const amrex::Real seed_y = static_cast<amrex::Real>(a_seed_y);
    const amrex::Real seed_z = static_cast<amrex::Real>(a_seed_z);

    const int c0 = a_fields.comps[0];
    const int c1 = a_fields.comps[1];
    const int c2 = a_fields.comps[2];
    const int c3 = a_fields.comps[3];

    // s = w^2 at (i,j,k), or -1 when the cell is outside the search ball.
    auto squared_amplitude =
        [=] AMREX_GPU_DEVICE(const amrex::Array4<const amrex::Real> &arr, int i,
                             int j, int k) -> amrex::Real
    {
        const amrex::Real px = prob_lo[0] + (amrex::Real(i) + 0.5) * dx_arr[0];
        const amrex::Real py = prob_lo[1] + (amrex::Real(j) + 0.5) * dx_arr[1];
        const amrex::Real pz = prob_lo[2] + (amrex::Real(k) + 0.5) * dx_arr[2];
        const amrex::Real ddx = px - seed_x;
        const amrex::Real ddy = py - seed_y;
        const amrex::Real ddz = pz - seed_z;
        if (ddx * ddx + ddy * ddy + ddz * ddz > ball_r2)
        {
            return amrex::Real(-1.0);
        }
        amrex::Real s = 0.0;
        if (c0 >= 0)
        {
            const amrex::Real v = arr(i, j, k, c0);
            s += v * v;
        }
        if (c1 >= 0)
        {
            const amrex::Real v = arr(i, j, k, c1);
            s += v * v;
        }
        if (c2 >= 0)
        {
            const amrex::Real v = arr(i, j, k, c2);
            s += v * v;
        }
        if (c3 >= 0)
        {
            const amrex::Real v = arr(i, j, k, c3);
            s += v * v;
        }
        return s;
    };

    // --- pass 1: the peak, which sets the core threshold -------------------
    amrex::Real peak_s = 0.0;
    {
        amrex::ReduceOps<amrex::ReduceOpMax> ops;
        amrex::ReduceData<amrex::Real> data(ops);
        using Tuple = typename decltype(data)::Type;
        for (amrex::MFIter mfi(a_state, amrex::TilingIfNotGPU()); mfi.isValid();
             ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            const auto arr       = a_state.const_array(mfi);
            ops.eval(bx, data,
                     [=] AMREX_GPU_DEVICE(int i, int j, int k) -> Tuple
                     {
                         const amrex::Real s = squared_amplitude(arr, i, j, k);
                         return {s > 0.0 ? s : amrex::Real(0.0)};
                     });
        }
        peak_s = amrex::get<0>(data.value());
        amrex::ParallelDescriptor::ReduceRealMax(peak_s);
    }

    SectorCore out;
    if (!(peak_s > 0.0))
    {
        return out; // no matter here: found stays false, caller holds its value
    }
    out.peak = std::sqrt(static_cast<double>(peak_s));

    // --- pass 2: the w^4-weighted centroid of the core ---------------------
    const amrex::Real core_s =
        static_cast<amrex::Real>(a_core_fraction * a_core_fraction) * peak_s;
    {
        amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpSum,
                         amrex::ReduceOpSum, amrex::ReduceOpSum>
            ops;
        amrex::ReduceData<amrex::Real, amrex::Real, amrex::Real, amrex::Real>
            data(ops);
        using Tuple = typename decltype(data)::Type;
        for (amrex::MFIter mfi(a_state, amrex::TilingIfNotGPU()); mfi.isValid();
             ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            const auto arr       = a_state.const_array(mfi);
            ops.eval(
                bx, data,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) -> Tuple
                {
                    const amrex::Real s = squared_amplitude(arr, i, j, k);
                    if (s < core_s)
                    {
                        return {0.0, 0.0, 0.0, 0.0};
                    }
                    const amrex::Real w4 = s * s;
                    const amrex::Real px =
                        prob_lo[0] + (amrex::Real(i) + 0.5) * dx_arr[0];
                    const amrex::Real py =
                        prob_lo[1] + (amrex::Real(j) + 0.5) * dx_arr[1];
                    const amrex::Real pz =
                        prob_lo[2] + (amrex::Real(k) + 0.5) * dx_arr[2];
                    return {w4, w4 * px, w4 * py, w4 * pz};
                });
        }
        auto v             = data.value();
        amrex::Real sum_w  = amrex::get<0>(v);
        amrex::Real sum_wx = amrex::get<1>(v);
        amrex::Real sum_wy = amrex::get<2>(v);
        amrex::Real sum_wz = amrex::get<3>(v);
        amrex::ParallelDescriptor::ReduceRealSum(sum_w);
        amrex::ParallelDescriptor::ReduceRealSum(sum_wx);
        amrex::ParallelDescriptor::ReduceRealSum(sum_wy);
        amrex::ParallelDescriptor::ReduceRealSum(sum_wz);

        if (sum_w > 0.0)
        {
            out.x     = static_cast<double>(sum_wx / sum_w);
            out.y     = static_cast<double>(sum_wy / sum_w);
            out.z     = static_cast<double>(sum_wz / sum_w);
            out.found = true;
        }
    }
    return out;
}

#endif /* SECTORCORETRACKER_HPP_ */
