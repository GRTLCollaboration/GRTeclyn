#ifndef RL_LUMP_TRACKER_HPP_
#define RL_LUMP_TRACKER_HPP_

#include "GRAMR.hpp"
#include "RLLumpState.hpp"
#include "RadialRecipeLevel.hpp"
#include "SetupFunctions.hpp"
#include "StateVariables.hpp"

#include <AMReX_Reduce.H>
#include <array>
#include <cmath>

//! Update per-lump 3-D state by localized reductions on the finest level.
//!
//! On entry, ``lumps[k].{x,y,z}`` hold each lump's previously tracked centre
//! (used as the search-ball centre); on exit every field of
//! ``lumps[0..num_lumps-1]`` is refreshed and the centre has moved to the new
//! mass-weighted centroid (so the matching pump "spotlight" follows the lump).
//!
//! Weight field = |Phi|^2.  For the complex boson field every lump is a blob of
//! the same field (distinguished only by its ball location); for independent
//! real scalars lump ``k`` uses its own ``phi_lump(k)`` component.
inline void track_rl_lumps(GRAMR &amr, int num_lumps, bool complex_field,
                           double ball_radius,
                           const std::array<double, AMREX_SPACEDIM> &center,
                           std::array<RLLumpState, RL_MAX_LUMPS> &lumps)
{
    if (num_lumps < 1)
        return;
    if (num_lumps > RL_MAX_LUMPS)
        num_lumps = RL_MAX_LUMPS;

    const int finest_lev = amr.finestLevel();
    auto &fine_level =
        dynamic_cast<RadialRecipeLevel &>(amr.getLevel(finest_lev));
    amrex::MultiFab &state_fine =
        fine_level.get_new_data(RadialRecipeLevel::state_index);
    const amrex::Real time =
        fine_level.get_state_data(RadialRecipeLevel::state_index).curTime();
    FillPatch(fine_level, state_fine, 2, time, RadialRecipeLevel::state_index, 0,
              state_fine.nComp());

    const auto geom = amr.Geom(finest_lev);
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx_arr =
        geom.CellSizeArray();
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_lo =
        geom.ProbLoArray();
    const amrex::Real cell_vol =
        AMREX_D_TERM(dx_arr[0], *dx_arr[1], *dx_arr[2]);
    const amrex::Real ball_r2 =
        static_cast<amrex::Real>(ball_radius * ball_radius);
    // Work in the same frame as the pump: physical coords relative to the grid
    // centre (Coordinates(..., grid_center)).  Lump centres / seeds are offsets.
    const amrex::Real cgx = static_cast<amrex::Real>(center[0]);
    const amrex::Real cgy = static_cast<amrex::Real>(center[1]);
    const amrex::Real cgz = static_cast<amrex::Real>(center[2]);

    for (int lump = 0; lump < num_lumps; ++lump)
    {
        const amrex::Real cx = static_cast<amrex::Real>(lumps[lump].x);
        const amrex::Real cy = static_cast<amrex::Real>(lumps[lump].y);
        const amrex::Real cz = static_cast<amrex::Real>(lumps[lump].z);
        const int phi_comp  = complex_field ? c_phi : c_phi_lump_index(lump);
        const int phi2_comp = complex_field ? c_phi2 : -1;

        amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpSum,
                         amrex::ReduceOpSum, amrex::ReduceOpSum,
                         amrex::ReduceOpSum, amrex::ReduceOpMax,
                         amrex::ReduceOpMin, amrex::ReduceOpMin>
            ops;
        amrex::ReduceData<amrex::Real, amrex::Real, amrex::Real, amrex::Real,
                          amrex::Real, amrex::Real, amrex::Real, amrex::Real>
            data(ops);
        using RT = typename decltype(data)::Type;

        for (amrex::MFIter mfi(state_fine, amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            const auto arr       = state_fine.const_array(mfi);
            ops.eval(
                bx, data,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) -> RT
                {
                    const amrex::Real px =
                        prob_lo[0] + (amrex::Real(i) + 0.5) * dx_arr[0] - cgx;
                    const amrex::Real py =
                        prob_lo[1] + (amrex::Real(j) + 0.5) * dx_arr[1] - cgy;
                    const amrex::Real pz =
                        prob_lo[2] + (amrex::Real(k) + 0.5) * dx_arr[2] - cgz;
                    const amrex::Real ddx = px - cx;
                    const amrex::Real ddy = py - cy;
                    const amrex::Real ddz = pz - cz;
                    const amrex::Real d2 = ddx * ddx + ddy * ddy + ddz * ddz;
                    if (d2 > ball_r2)
                    {
                        // Out of ball: neutral for every reduction op.
                        return {0.0,    0.0,    0.0,   0.0,
                                0.0,    0.0,    1.0e30, 1.0e30};
                    }
                    amrex::Real w = arr(i, j, k, phi_comp);
                    w *= w;
                    if (phi2_comp >= 0)
                    {
                        const amrex::Real p2 = arr(i, j, k, phi2_comp);
                        w += p2 * p2;
                    }
                    const amrex::Real lapse = arr(i, j, k, c_lapse);
                    const amrex::Real chi   = arr(i, j, k, c_chi);
                    return {w,   w * ddx, w * ddy, w * ddz,
                            w * d2, w,    lapse,   chi};
                });
        }

        auto hv             = data.value();
        amrex::Real sum_w   = amrex::get<0>(hv);
        amrex::Real sum_wdx = amrex::get<1>(hv);
        amrex::Real sum_wdy = amrex::get<2>(hv);
        amrex::Real sum_wdz = amrex::get<3>(hv);
        amrex::Real sum_wr2 = amrex::get<4>(hv);
        amrex::Real max_w   = amrex::get<5>(hv);
        amrex::Real min_lap = amrex::get<6>(hv);
        amrex::Real min_chi = amrex::get<7>(hv);
        amrex::ParallelDescriptor::ReduceRealSum(sum_w);
        amrex::ParallelDescriptor::ReduceRealSum(sum_wdx);
        amrex::ParallelDescriptor::ReduceRealSum(sum_wdy);
        amrex::ParallelDescriptor::ReduceRealSum(sum_wdz);
        amrex::ParallelDescriptor::ReduceRealSum(sum_wr2);
        amrex::ParallelDescriptor::ReduceRealMax(max_w);
        amrex::ParallelDescriptor::ReduceRealMin(min_lap);
        amrex::ParallelDescriptor::ReduceRealMin(min_chi);

        RLLumpState &out = lumps[lump];
        out.min_lapse =
            (min_lap < 1.0e29) ? static_cast<double>(min_lap) : 1.0;
        out.min_chi = (min_chi < 1.0e29) ? static_cast<double>(min_chi) : 1.0;
        if (sum_w > 0.0)
        {
            const double ox = sum_wdx / sum_w; // centroid offset from prev centre
            const double oy = sum_wdy / sum_w;
            const double oz = sum_wdz / sum_w;
            out.x           = static_cast<double>(cx) + ox;
            out.y           = static_cast<double>(cy) + oy;
            out.z           = static_cast<double>(cz) + oz;
            const double var =
                sum_wr2 / sum_w - (ox * ox + oy * oy + oz * oz);
            out.size = (var > 0.0) ? std::sqrt(var) : 0.0;
            out.mass = static_cast<double>(sum_w * cell_vol);
            out.peak = static_cast<double>(max_w);
        }
        else
        {
            // Lump dispersed below detectability: hold last centre, flag empty.
            out.size = ball_radius;
            out.mass = 0.0;
            out.peak = 0.0;
        }
    }
}

#endif /* RL_LUMP_TRACKER_HPP_ */
