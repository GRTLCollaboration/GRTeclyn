/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef GRIDTREADMILL_HPP_
#define GRIDTREADMILL_HPP_

#include <AMReX_Geometry.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MultiFab.H>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

/// Exact whole-cell translation of the evolved state ("recentring box").
///
/// A self-accelerating source eventually walks out of a fixed box.  Rather than
/// enlarging the box -- whose cost grows with the target speed -- this module
/// periodically carries the *data* back toward the box centre by a whole number
/// of cells and accumulates the distance travelled in an odometer:
///
///     true displacement = position on the grid + odometer
///
/// Because the translation is a whole number of cells it is a pure relabelling
/// of which cell holds which value: no interpolation, exact to the last bit.
/// A translation that landed between cells would have to average neighbouring
/// values and would quietly bias the run -- so only whole-cell shifts are ever
/// taken, and `cells_for_length` refuses anything else.
///
/// One seam is unavoidable, at the face the data is carried away from:
///
///   * that face loses `|n|` layers of source and must be refilled;
///   * the opposite face has `|n|` layers of already-departed wake pushed out
///     of the domain and discarded.
///
/// Both live in the outermost shell, where the sponge (`SpongeZone.hpp`)
/// absorbs whatever is there.  The caller is responsible for checking that the
/// shift length keeps the source well inside the sponge's inner radius.
///
/// SINGLE RESPONSIBILITY.  This file knows nothing about any particular
/// example: it takes a MultiFab, an axis, a signed cell count and a table of
/// asymptotic values.  The decision of *when* to shift, and what counts as the
/// source, belongs to the caller.
struct GridTreadmillParams
{
    bool enabled{false};
    int axis{0};             //!< 0/1/2 -- the axis the source runs along
    double threshold{2.0};   //!< shift once the source is this far from centre
    int check_interval{20};  //!< steps between position reductions
    double ball_radius{8.0}; //!< search-ball radius for the source tracker
    int fill_mode{0};        //!< 0 = asymptotic 1/r falloff, 1 = copy outermost
};

namespace GridTreadmill
{

/// Fill strategies for the invented sliver.
///
/// OPEN/CLOSED: a new strategy is a new value handled in `fill_sliver`; nothing
/// else in the module changes.  Having two that disagree is deliberate -- a
/// trajectory that depends on which is used is telling you the invented region
/// matters, which is exactly the thing that has to be tested rather than
/// assumed.
enum FillMode
{
    FILL_ASYMPTOTIC = 0, //!< v = v_inf + (v_src - v_inf) * r_src / r
    FILL_COPY       = 1  //!< v = v_src (zeroth-order extrapolation)
};

/// Whole number of cells corresponding to `a_length`, or -1 if `a_length` is
/// not a whole number of cells.
///
/// The -1 is deliberately not a fallback to something reasonable: a shift that
/// does not land on the mesh cannot be taken exactly, and silently rounding it
/// would reintroduce precisely the sub-cell aliasing this design exists to
/// avoid -- the same failure that displaced the solved metric relative to the
/// matter in the 2026-08-21 initial-data bug.
inline int cells_for_length(double a_length, double a_dx, double a_tol = 1.0e-9)
{
    if (!(a_dx > 0.0) || !(a_length > 0.0))
    {
        return -1;
    }
    const double cells   = a_length / a_dx;
    const double rounded = std::round(cells);
    if (rounded < 1.0 ||
        std::abs(cells - rounded) > a_tol * std::max(1.0, cells))
    {
        return -1;
    }
    return static_cast<int>(rounded);
}

/// Smallest extent, along `a_axis`, of any box in `a_ba`.
///
/// The sliver is filled from the outermost surviving layer, which is only in
/// the same box as the sliver when the shift is shorter than a box.  The caller
/// checks this once at startup rather than discovering it as a wrong answer on
/// one rank.
inline int min_box_extent(const amrex::BoxArray &a_ba, int a_axis)
{
    if (a_ba.size() == 0)
    {
        return 0;
    }
    int smallest = std::numeric_limits<int>::max();
    for (int i = 0; i < static_cast<int>(a_ba.size()); ++i)
    {
        smallest = std::min(smallest, a_ba[i].length(a_axis));
    }
    return smallest;
}

/// The index range that `translate` leaves without a source, for a signed shift.
///
/// Positive `a_cells` carries the data toward the low end of the axis, so the
/// high face is left empty; negative does the reverse.  Returned as a Box so
/// callers and tests refer to one definition of "the sliver" rather than each
/// deriving it.
inline amrex::Box sliver_box(const amrex::Box &a_domain, int a_axis,
                             int a_cells)
{
    amrex::Box out = a_domain;
    if (a_cells > 0)
    {
        out.setSmall(a_axis, a_domain.bigEnd(a_axis) - a_cells + 1);
    }
    else if (a_cells < 0)
    {
        out.setBig(a_axis, a_domain.smallEnd(a_axis) - a_cells - 1);
    }
    return out;
}

/// Index of the outermost surviving layer that the sliver is filled from.
inline int sliver_source_index(const amrex::Box &a_domain, int a_axis,
                               int a_cells)
{
    return (a_cells > 0) ? (a_domain.bigEnd(a_axis) - a_cells)
                         : (a_domain.smallEnd(a_axis) - a_cells);
}

/// Translate every cell of `a_state` by `a_cells` along `a_axis`, so that
/// afterwards
///
///     a_state_new(i) = a_state_old(i + a_cells)
///
/// on the overlap.  The `|a_cells|` layers named by `sliver_box` are left with
/// no source and MUST be refilled by `fill_sliver` -- this function does not do
/// it, so that "translate" and "invent the missing shell" stay separate
/// operations that can be tested and reasoned about one at a time.
///
/// The overlap is copied bit-for-bit.  Ghost cells are NOT filled -- the caller
/// re-fills them (any `FillPatch` does) before the next stencil evaluation.
inline void translate(amrex::MultiFab &a_state, int a_axis, int a_cells)
{
    AMREX_ALWAYS_ASSERT(a_cells != 0);
    AMREX_ALWAYS_ASSERT(a_axis >= 0 && a_axis < AMREX_SPACEDIM);

    const int ncomp = a_state.nComp();

    // FabArray::shift relabels the BoxArray and every FAB's index space
    // together, so no data moves within a FAB; ParallelCopy then delivers each
    // cell to its new owner.  Both steps are exact.
    amrex::MultiFab tmp(a_state.boxArray(), a_state.DistributionMap(), ncomp, 0);
    amrex::MultiFab::Copy(tmp, a_state, 0, 0, ncomp, 0);

    amrex::IntVect shift_iv(0);
    shift_iv[a_axis] = -a_cells;
    tmp.shift(shift_iv);

    // ParallelCopy writes ONLY the intersection.  The sliver layers have no
    // source and are left holding their pre-shift values -- stale data that
    // looks entirely plausible.  That is why filling them is a separate,
    // mandatory step and not an optional tidy-up.
    a_state.ParallelCopy(tmp, 0, 0, ncomp, 0, 0,
                         amrex::Periodicity::NonPeriodic());
}

/// Invent the layers that `translate` could not supply, from the outermost
/// surviving layer.
inline void fill_sliver(amrex::MultiFab &a_state, int a_axis, int a_cells,
                        const amrex::Geometry &a_geom,
                        const std::array<double, AMREX_SPACEDIM> &a_center,
                        const std::vector<double> &a_asymptotic_values,
                        int a_fill_mode)
{
    AMREX_ALWAYS_ASSERT(a_cells != 0);
    AMREX_ALWAYS_ASSERT(a_axis >= 0 && a_axis < AMREX_SPACEDIM);

    const int ncomp = a_state.nComp();
    AMREX_ALWAYS_ASSERT(static_cast<int>(a_asymptotic_values.size()) >= ncomp);

    const amrex::Box domain = a_geom.Domain();
    const amrex::Box sliver = sliver_box(domain, a_axis, a_cells);
    const int ref_index     = sliver_source_index(domain, a_axis, a_cells);

    amrex::Gpu::DeviceVector<amrex::Real> d_asymptotic(ncomp);
    {
        std::vector<amrex::Real> h_asymptotic(ncomp);
        for (int n = 0; n < ncomp; ++n)
        {
            h_asymptotic[n] = static_cast<amrex::Real>(a_asymptotic_values[n]);
        }
        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_asymptotic.begin(),
                              h_asymptotic.end(), d_asymptotic.begin());
        amrex::Gpu::streamSynchronize();
    }
    const amrex::Real *asymptotic = d_asymptotic.data();

    const auto prob_lo = a_geom.ProbLoArray();
    const auto dx_arr  = a_geom.CellSizeArray();
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> centre{
        AMREX_D_DECL(static_cast<amrex::Real>(a_center[0]),
                     static_cast<amrex::Real>(a_center[1]),
                     static_cast<amrex::Real>(a_center[2]))};
    const int axis       = a_axis;
    const int ref        = ref_index;
    const bool copy_fill = (a_fill_mode == FILL_COPY);

    for (amrex::MFIter mfi(a_state, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi)
    {
        const amrex::Box bx = mfi.validbox() & sliver;
        if (!bx.ok())
        {
            continue;
        }
        const auto arr = a_state.array(mfi);
        amrex::ParallelFor(
            bx, ncomp,
            [=] AMREX_GPU_DEVICE(int i, int j, int k, int n)
            {
                int si = i;
                int sj = j;
                int sk = k;
                if (axis == 0)
                {
                    si = ref;
                }
                else if (axis == 1)
                {
                    sj = ref;
                }
                else
                {
                    sk = ref;
                }

                const amrex::Real v_src = arr(si, sj, sk, n);
                if (copy_fill)
                {
                    arr(i, j, k, n) = v_src;
                    return;
                }

                const amrex::Real x =
                    prob_lo[0] + (amrex::Real(i) + 0.5) * dx_arr[0] - centre[0];
                const amrex::Real y =
                    prob_lo[1] + (amrex::Real(j) + 0.5) * dx_arr[1] - centre[1];
                const amrex::Real z =
                    prob_lo[2] + (amrex::Real(k) + 0.5) * dx_arr[2] - centre[2];
                const amrex::Real sx = prob_lo[0] +
                                       (amrex::Real(si) + 0.5) * dx_arr[0] -
                                       centre[0];
                const amrex::Real sy = prob_lo[1] +
                                       (amrex::Real(sj) + 0.5) * dx_arr[1] -
                                       centre[1];
                const amrex::Real sz = prob_lo[2] +
                                       (amrex::Real(sk) + 0.5) * dx_arr[2] -
                                       centre[2];

                const amrex::Real r = std::sqrt(x * x + y * y + z * z);
                const amrex::Real r_src =
                    std::sqrt(sx * sx + sy * sy + sz * sz);

                const amrex::Real v_inf = asymptotic[n];
                if (r > 0.0)
                {
                    // The same 1/r approach to the asymptotic value that the
                    // Sommerfeld boundary condition already assumes, so the
                    // fabricated cells are continuous with what the boundary
                    // would have produced anyway.
                    arr(i, j, k, n) = v_inf + (v_src - v_inf) * (r_src / r);
                }
                else
                {
                    arr(i, j, k, n) = v_src;
                }
            });
    }
    amrex::Gpu::streamSynchronize();
}

/// Translate and refill, in the only order that is ever correct.
///
/// Provided so no caller can accidentally do one without the other; the two
/// steps remain separately callable for the unit tests, which have to look at
/// the state in between.
inline void shift(amrex::MultiFab &a_state, int a_axis, int a_cells,
                  const amrex::Geometry &a_geom,
                  const std::array<double, AMREX_SPACEDIM> &a_center,
                  const std::vector<double> &a_asymptotic_values,
                  int a_fill_mode)
{
    translate(a_state, a_axis, a_cells);
    fill_sliver(a_state, a_axis, a_cells, a_geom, a_center, a_asymptotic_values,
                a_fill_mode);
}

/// The odometer: how far the data has been carried, in cells and in length.
///
/// This is the only genuinely new piece of run state, and it is the one thing
/// that MUST survive a restart.  Without it a restarted run resumes with the
/// trajectory silently reset to zero and produces a plot that looks completely
/// plausible while being wrong.
struct Odometer
{
    long cells{0};
    double dx{0.0};

    [[nodiscard]] double length() const
    {
        return static_cast<double>(cells) * dx;
    }

    void add(int a_cells) { cells += a_cells; }

    /// Serialised form written beside a checkpoint.  Plain text on purpose: it
    /// has to be readable by a human deciding whether a restart is sane.
    [[nodiscard]] std::string to_string() const
    {
        std::ostringstream os;
        os.precision(17);
        os << "treadmill_odometer_cells " << cells << "\n"
           << "treadmill_dx " << dx << "\n";
        return os.str();
    }
};

} // namespace GridTreadmill

#endif /* GRIDTREADMILL_HPP_ */
