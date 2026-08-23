/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// Doctest header
#include "doctest.h"

// Test header
#include "GridTreadmillTest.hpp"

// Common includes
#include "doctestCLIArgs.hpp"

// AMReX includes
#include "AMReX.H"
#include "AMReX_MultiFab.H"

// Other includes
#include "GridTreadmill.hpp"

#include <cmath>
#include <sstream>
#include <vector>

namespace
{

constexpr int num_cells = 16;
constexpr int num_comps = 3;
constexpr int max_box   = 8;
constexpr double dx     = 1.0;

//! A value that is unique per (cell, component), so a mis-addressed copy cannot
//! coincidentally look right.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real pattern(int i, int j,
                                                             int k, int n)
{
    return static_cast<amrex::Real>(i) + 100.0 * static_cast<amrex::Real>(j) +
           10000.0 * static_cast<amrex::Real>(k) +
           1000000.0 * static_cast<amrex::Real>(n);
}

struct TestGrid
{
    amrex::Box domain;
    amrex::RealBox real_box;
    amrex::Geometry geom;
    amrex::BoxArray box_array;
    amrex::DistributionMapping dm;
    amrex::MFInfo mf_info;

    TestGrid()
        : domain{amrex::IntVect::TheZeroVector(),
                 amrex::IntVect{num_cells - 1}},
          real_box{}, geom{}, box_array{}, dm{}, mf_info{}
    {
        amrex::RealVect dx_vect{dx};
        real_box = amrex::RealBox{domain, dx_vect.dataPtr(),
                                  amrex::RealVect::Zero.dataPtr()};
        geom     = amrex::Geometry{domain, &real_box, 0};
        box_array = amrex::BoxArray{domain};
        box_array.maxSize(max_box);
        dm = amrex::DistributionMapping{box_array};
        mf_info.SetArena(amrex::The_Managed_Arena());
    }

    [[nodiscard]] amrex::MultiFab make_patterned() const
    {
        amrex::MultiFab mf{box_array, dm, num_comps, 0, mf_info};
        const auto &arrays = mf.arrays();
        amrex::ParallelFor(mf, amrex::IntVect(0), num_comps,
                           [=] AMREX_GPU_DEVICE(int ibox, int i, int j, int k,
                                                int n)
                           { arrays[ibox](i, j, k, n) = pattern(i, j, k, n); });
        amrex::Gpu::streamSynchronize();
        return mf;
    }
};

//! Asymptotic values: distinct per component so the fill cannot accidentally
//! agree with the pattern.
std::vector<double> asymptotic_values()
{
    return {0.5, -1.25, 3.0};
}

std::array<double, AMREX_SPACEDIM> box_centre()
{
    const double half = 0.5 * num_cells * dx;
    return {AMREX_D_DECL(half, half, half)};
}

} // namespace

void run_grid_treadmill_test()
{
    int amrex_argc    = doctest::cli_args.argc();
    char **amrex_argv = doctest::cli_args.argv();
    // NOLINTNEXTLINE(bugprone-casting-through-void) // Open MPI triggers this
    amrex::Initialize(amrex_argc, amrex_argv);
    {
        const TestGrid grid;
        const auto centre     = box_centre();
        const auto asymptotic = asymptotic_values();

        // -- a shift must be a whole number of cells, or it is refused -------
        //
        // Rounding a sub-cell shift would reintroduce exactly the aliasing this
        // module exists to avoid, so "not representable" has to be an error and
        // not a nearest-neighbour guess.
        SUBCASE("cells_for_length refuses anything that is not whole cells")
        {
            CHECK(GridTreadmill::cells_for_length(2.0, 0.5) == 4);
            CHECK(GridTreadmill::cells_for_length(0.5, 0.5) == 1);
            CHECK(GridTreadmill::cells_for_length(0.75, 0.5) == -1);
            CHECK(GridTreadmill::cells_for_length(2.0, 0.3) == -1);
            CHECK(GridTreadmill::cells_for_length(0.0, 0.5) == -1);
            CHECK(GridTreadmill::cells_for_length(-1.0, 0.5) == -1);
        }

        // -- the translation is exact on the overlap -------------------------
        SUBCASE("translate reproduces the pattern bit-for-bit")
        {
            constexpr int shift = 2;
            amrex::MultiFab mf  = grid.make_patterned();
            GridTreadmill::translate(mf, 0, shift);
            amrex::Gpu::streamSynchronize();

            const amrex::Box overlap = amrex::Box{
                grid.domain.smallEnd(),
                amrex::IntVect{grid.domain.bigEnd(0) - shift,
                               grid.domain.bigEnd(1), grid.domain.bigEnd(2)}};

            int mismatches = 0;
            for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi)
            {
                const amrex::Box bx = mfi.validbox() & overlap;
                if (!bx.ok())
                {
                    continue;
                }
                const auto arr = mf.const_array(mfi);
                amrex::LoopOnCpu(bx, num_comps,
                                 [&](int i, int j, int k, int n)
                                 {
                                     // Bit-for-bit: a tolerance here would hide
                                     // precisely the interpolation this design
                                     // promises never to do.
                                     if (arr(i, j, k, n) !=
                                         pattern(i + shift, j, k, n))
                                     {
                                         ++mismatches;
                                     }
                                 });
            }
            CHECK(mismatches == 0);
        }

        // -- the sliver is always overwritten --------------------------------
        //
        // ParallelCopy writes only the intersection, so without an explicit
        // fill the leading layers keep their pre-shift values: a plausible
        // trajectory built on stale data.  This is the likeliest bug in the
        // module and the reason the test exists.
        SUBCASE("fill_sliver leaves no pre-shift value behind")
        {
            constexpr int shift          = 2;
            constexpr amrex::Real poison = -987654.0;
            amrex::MultiFab mf           = grid.make_patterned();

            const amrex::Box sliver =
                GridTreadmill::sliver_box(grid.domain, 0, shift);
            {
                const auto &arrays = mf.arrays();
                amrex::ParallelFor(
                    mf, amrex::IntVect(0), num_comps,
                    [=] AMREX_GPU_DEVICE(int ibox, int i, int j, int k, int n)
                    {
                        if (sliver.contains(amrex::IntVect{i, j, k}))
                        {
                            arrays[ibox](i, j, k, n) = poison;
                        }
                    });
                amrex::Gpu::streamSynchronize();
            }

            GridTreadmill::translate(mf, 0, shift);
            GridTreadmill::fill_sliver(mf, 0, shift, grid.geom, centre,
                                       asymptotic,
                                       GridTreadmill::FILL_ASYMPTOTIC);
            amrex::Gpu::streamSynchronize();

            int survivors = 0;
            for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi)
            {
                const amrex::Box bx = mfi.validbox() & sliver;
                if (!bx.ok())
                {
                    continue;
                }
                const auto arr = mf.const_array(mfi);
                amrex::LoopOnCpu(bx, num_comps,
                                 [&](int i, int j, int k, int n)
                                 {
                                     if (arr(i, j, k, n) == poison)
                                     {
                                         ++survivors;
                                     }
                                 });
            }
            CHECK(survivors == 0);
        }

        // -- both fill strategies stay between the source and the asymptote --
        SUBCASE("the asymptotic fill relaxes toward the asymptotic value")
        {
            constexpr int shift = 2;
            amrex::MultiFab copy_filled = grid.make_patterned();
            amrex::MultiFab asym_filled = grid.make_patterned();

            GridTreadmill::shift(copy_filled, 0, shift, grid.geom, centre,
                                 asymptotic, GridTreadmill::FILL_COPY);
            GridTreadmill::shift(asym_filled, 0, shift, grid.geom, centre,
                                 asymptotic, GridTreadmill::FILL_ASYMPTOTIC);
            amrex::Gpu::streamSynchronize();

            const amrex::Box sliver =
                GridTreadmill::sliver_box(grid.domain, 0, shift);

            int out_of_range = 0;
            for (amrex::MFIter mfi(asym_filled); mfi.isValid(); ++mfi)
            {
                const amrex::Box bx = mfi.validbox() & sliver;
                if (!bx.ok())
                {
                    continue;
                }
                const auto a = asym_filled.const_array(mfi);
                const auto c = copy_filled.const_array(mfi);
                amrex::LoopOnCpu(
                    bx, num_comps,
                    [&](int i, int j, int k, int n)
                    {
                        // The 1/r fill moves the copied value toward the
                        // asymptote and never past it, so it must lie in the
                        // closed interval between the two.
                        const double lo =
                            std::min(static_cast<double>(c(i, j, k, n)),
                                     asymptotic[n]);
                        const double hi =
                            std::max(static_cast<double>(c(i, j, k, n)),
                                     asymptotic[n]);
                        const double v = a(i, j, k, n);
                        if (v < lo - 1.0e-9 || v > hi + 1.0e-9)
                        {
                            ++out_of_range;
                        }
                    });
            }
            CHECK(out_of_range == 0);
        }

        // -- forward then back restores the interior exactly ------------------
        //
        // This is the unit-level version of the forced-shift null cell: the
        // seams are exercised twice and the physics region must come back
        // untouched.
        SUBCASE("a there-and-back shift leaves the interior bit-identical")
        {
            constexpr int shift = 2;
            amrex::MultiFab mf  = grid.make_patterned();

            GridTreadmill::shift(mf, 0, shift, grid.geom, centre, asymptotic,
                                 GridTreadmill::FILL_ASYMPTOTIC);
            GridTreadmill::shift(mf, 0, -shift, grid.geom, centre, asymptotic,
                                 GridTreadmill::FILL_ASYMPTOTIC);
            amrex::Gpu::streamSynchronize();

            // Both seams are touched, so only the region that neither shift
            // invented is required to be untouched.
            const amrex::Box interior = amrex::Box{
                amrex::IntVect{grid.domain.smallEnd(0) + shift,
                               grid.domain.smallEnd(1), grid.domain.smallEnd(2)},
                amrex::IntVect{grid.domain.bigEnd(0) - shift,
                               grid.domain.bigEnd(1), grid.domain.bigEnd(2)}};

            int mismatches = 0;
            for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi)
            {
                const amrex::Box bx = mfi.validbox() & interior;
                if (!bx.ok())
                {
                    continue;
                }
                const auto arr = mf.const_array(mfi);
                amrex::LoopOnCpu(bx, num_comps,
                                 [&](int i, int j, int k, int n)
                                 {
                                     if (arr(i, j, k, n) != pattern(i, j, k, n))
                                     {
                                         ++mismatches;
                                     }
                                 });
            }
            CHECK(mismatches == 0);
        }

        // -- a shift toward the low face is the mirror of one toward the high -
        SUBCASE("a negative shift empties the low face instead")
        {
            constexpr int shift = 3;
            amrex::MultiFab mf  = grid.make_patterned();
            GridTreadmill::translate(mf, 0, -shift);
            amrex::Gpu::streamSynchronize();

            const amrex::Box overlap = amrex::Box{
                amrex::IntVect{grid.domain.smallEnd(0) + shift,
                               grid.domain.smallEnd(1),
                               grid.domain.smallEnd(2)},
                grid.domain.bigEnd()};

            int mismatches = 0;
            for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi)
            {
                const amrex::Box bx = mfi.validbox() & overlap;
                if (!bx.ok())
                {
                    continue;
                }
                const auto arr = mf.const_array(mfi);
                amrex::LoopOnCpu(bx, num_comps,
                                 [&](int i, int j, int k, int n)
                                 {
                                     if (arr(i, j, k, n) !=
                                         pattern(i - shift, j, k, n))
                                     {
                                         ++mismatches;
                                     }
                                 });
            }
            CHECK(mismatches == 0);
            CHECK(GridTreadmill::sliver_box(grid.domain, 0, -shift).bigEnd(0) ==
                  grid.domain.smallEnd(0) + shift - 1);
        }

        // -- the odometer survives being written down and read back ----------
        //
        // A restart that resumes with the odometer at zero produces a
        // trajectory that is wrong and looks entirely plausible, so the
        // serialised form is checked rather than assumed.
        SUBCASE("the odometer round-trips through its serialised form")
        {
            GridTreadmill::Odometer odo;
            odo.dx = 0.5;
            odo.add(4);
            odo.add(4);
            odo.add(-4);
            CHECK(odo.cells == 4);
            CHECK(odo.length() == doctest::Approx(2.0));

            std::istringstream is(odo.to_string());
            std::string key;
            long cells = -1;
            double read_dx = -1.0;
            while (is >> key)
            {
                if (key == "treadmill_odometer_cells")
                {
                    is >> cells;
                }
                else if (key == "treadmill_dx")
                {
                    is >> read_dx;
                }
            }
            CHECK(cells == odo.cells);
            CHECK(read_dx == doctest::Approx(odo.dx));
        }

        // -- the startup guard sees the real box decomposition ---------------
        SUBCASE("min_box_extent reports the smallest box on the axis")
        {
            CHECK(GridTreadmill::min_box_extent(grid.box_array, 0) == max_box);
        }
    }
    amrex::Finalize();
}
