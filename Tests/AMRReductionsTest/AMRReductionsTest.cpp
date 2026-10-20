/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// Doctest header
#include "doctest.h"

// Test header
#include "AMRReductionsTest.hpp"

// Common Test headers
#include "doctestCLIArgs.hpp"

// Our includes
#include "AMRReductions.hpp"

// AMReX headers
#include <AMReX.H>
#include <AMReX_MultiFab.H>

// System headers
#include <cmath>
#include <memory>

namespace
{
// Components of the test data
enum
{
    c_x,        // f = x
    c_constant, // f = -2
    c_vector_0, // (3, 4) has magnitude 5
    c_vector_1,
    NUM_TEST_VARS
};

constexpr amrex::Real constant_value = -2.0;

/* Builds a two-level hierarchy on the unit cube with num_coarse_cells^3
 * coarse cells. The fine level (refinement ratio 2) covers the corner
 * [0.5, 1]^3. Ghost cells are filled with a huge value that would spoil
 * any reduction they leak into. */
struct TestHierarchy
{
    static constexpr int num_coarse_cells = 16;
    static constexpr int max_grid_size    = 8;
    static constexpr int num_ghosts       = 1;

    amrex::Vector<amrex::Geometry> geometry;
    amrex::Vector<amrex::IntVect> refinement_ratios{amrex::IntVect(2)};
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> data;

    TestHierarchy()
    {
        const amrex::RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                                      {AMREX_D_DECL(1.0, 1.0, 1.0)});
        const amrex::Array<int, AMREX_SPACEDIM> is_periodic{
            AMREX_D_DECL(0, 0, 0)};
        const amrex::Box coarse_domain(amrex::IntVect(0),
                                       amrex::IntVect(num_coarse_cells - 1));
        const amrex::Box fine_domain = amrex::refine(coarse_domain, 2);
        geometry.emplace_back(coarse_domain, real_box, 0, is_periodic);
        geometry.emplace_back(fine_domain, real_box, 0, is_periodic);

        const amrex::Box fine_region(amrex::IntVect(num_coarse_cells),
                                     amrex::IntVect(2 * num_coarse_cells - 1));
        for (const amrex::Box &box : {coarse_domain, fine_region})
        {
            amrex::BoxArray box_array(box);
            box_array.maxSize(max_grid_size);
            const amrex::DistributionMapping distribution_mapping(box_array);
            data.push_back(std::make_unique<amrex::MultiFab>(
                box_array, distribution_mapping, NUM_TEST_VARS, num_ghosts));
        }

        for (int level = 0; level < data.size(); ++level)
        {
            auto &level_data = *data[level];
            level_data.setVal(1e30);
            const auto dx      = geometry[level].CellSizeArray();
            const auto &arrays = level_data.arrays();
            amrex::ParallelFor(
                level_data,
                // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy,
                                     int iz) noexcept
                {
                    const auto &cell             = arrays[box_no];
                    cell(ix, iy, iz, c_x)        = (ix + 0.5) * dx[0];
                    cell(ix, iy, iz, c_constant) = constant_value;
                    cell(ix, iy, iz, c_vector_0) = 3.0;
                    cell(ix, iy, iz, c_vector_1) = 4.0;
                });
        }
        amrex::Gpu::streamSynchronize();
    }
};
} // namespace

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
void run_amr_reductions_test()
{
    using AMRReductions::Integral;
    using AMRReductions::L1Norm;
    using AMRReductions::L2Norm;
    using AMRReductions::Max;
    using AMRReductions::Min;
    using AMRReductions::PNorm;

    int amrex_argc    = doctest::cli_args.argc();
    char **amrex_argv = doctest::cli_args.argv();
    // NOLINTNEXTLINE(bugprone-casting-through-void) // Open MPI triggers this
    amrex::Initialize(amrex_argc, amrex_argv);
    {
        const TestHierarchy hierarchy;
        const auto data = amrex::GetVecOfConstPtrs(hierarchy.data);
        constexpr amrex::Real fine_dx   = 0.5 / TestHierarchy::num_coarse_cells;
        constexpr amrex::Real coarse_dx = 1.0 / TestHierarchy::num_coarse_cells;
        const amrex::Real abs_constant  = std::abs(constant_value);

        SUBCASE("Every operation in a single call")
        {
            const AMRReductions::AMRReducer reducer(
                hierarchy.geometry, hierarchy.refinement_ratios);
            CHECK(reducer.domain_volume() == doctest::Approx(1.0));

            const auto [x_integral, x_min, x_max, constant_integral,
                        constant_l2, constant_l1, constant_l3, constant_min,
                        vector_l2, vector_integral, vector_max] =
                reducer.compute(
                    data, Integral{.comp = c_x}, Min{.comp = c_x},
                    Max{.comp = c_x}, Integral{.comp = c_constant},
                    L2Norm{.comp = c_constant}, L1Norm{.comp = c_constant},
                    PNorm<3>{.comp = c_constant}, Min{.comp = c_constant},
                    L2Norm{.comp = c_vector_0, .num_comps = 2},
                    Integral{.comp = c_vector_0, .num_comps = 2},
                    Max{.comp = c_vector_0, .num_comps = 2});

            // The midpoint rule is exact for a linear function, so this
            // would fail if covered coarse cells were counted twice
            CHECK(x_integral == doctest::Approx(0.5));
            // The minimum lies on the coarse level, the maximum on the fine
            CHECK(x_min == doctest::Approx(0.5 * coarse_dx));
            CHECK(x_max == doctest::Approx(1.0 - 0.5 * fine_dx));
            CHECK(constant_integral == doctest::Approx(constant_value));
            CHECK(constant_l2 == doctest::Approx(abs_constant));
            CHECK(constant_l1 == doctest::Approx(abs_constant));
            CHECK(constant_l3 == doctest::Approx(abs_constant));
            CHECK(constant_min == doctest::Approx(constant_value));
            CHECK(vector_l2 == doctest::Approx(5.0));
            CHECK(vector_integral == doctest::Approx(5.0));
            CHECK(vector_max == doctest::Approx(5.0));
        }

        SUBCASE("Unnormalized norm")
        {
            const AMRReductions::AMRReducer reducer(
                hierarchy.geometry, hierarchy.refinement_ratios);
            const auto [x_l1] = reducer.compute(
                data, L1Norm{.comp = c_x, .normalize_by_volume = false});
            CHECK(x_l1 == doctest::Approx(0.5));
        }

        SUBCASE("Excluding boundary cells")
        {
            // 2 coarse cells (i.e. 4 fine cells) are removed at every face,
            // leaving the cube [0.125, 0.875]^3
            constexpr int cells_to_exclude = 2;
            const AMRReductions::AMRReducer reducer(
                hierarchy.geometry, hierarchy.refinement_ratios,
                {.lo = amrex::IntVect(cells_to_exclude),
                 .hi = amrex::IntVect(cells_to_exclude)});
            constexpr amrex::Real included_lo = cells_to_exclude * coarse_dx;
            constexpr amrex::Real included_hi = 1.0 - included_lo;
            const amrex::Real volume =
                std::pow(included_hi - included_lo, AMREX_SPACEDIM);
            CHECK(reducer.domain_volume() == doctest::Approx(volume));

            const auto [x_integral, x_min, x_max, constant_l2] =
                reducer.compute(data, Integral{.comp = c_x}, Min{.comp = c_x},
                                Max{.comp = c_x}, L2Norm{.comp = c_constant});
            CHECK(x_integral == doctest::Approx(0.5 * volume));
            CHECK(x_min == doctest::Approx(included_lo + 0.5 * coarse_dx));
            CHECK(x_max == doctest::Approx(included_hi - 0.5 * fine_dx));
            CHECK(constant_l2 == doctest::Approx(abs_constant));
        }

        SUBCASE("Excluding boundary cells at a single face")
        {
            // Only the high x face loses 2 coarse cells (i.e. 4 fine cells),
            // leaving [0, 0.875] x [0, 1] x [0, 1]
            constexpr int cells_to_exclude = 2;
            const AMRReductions::AMRReducer reducer(
                hierarchy.geometry, hierarchy.refinement_ratios,
                {.hi = amrex::IntVect(AMREX_D_DECL(cells_to_exclude, 0, 0))});
            constexpr amrex::Real x_hi = 1.0 - cells_to_exclude * coarse_dx;
            CHECK(reducer.domain_volume() == doctest::Approx(x_hi));

            const auto [x_integral, x_min, x_max] =
                reducer.compute(data, Integral{.comp = c_x}, Min{.comp = c_x},
                                Max{.comp = c_x});
            CHECK(x_integral == doctest::Approx(0.5 * x_hi * x_hi));
            CHECK(x_min == doctest::Approx(0.5 * coarse_dx));
            CHECK(x_max == doctest::Approx(x_hi - 0.5 * fine_dx));
        }

        SUBCASE("Reductions in one call agree with separate calls")
        {
            const AMRReductions::AMRReducer reducer(
                hierarchy.geometry, hierarchy.refinement_ratios);
            const PNorm<3> l3_norm{.comp = c_x};
            const auto together = reducer.compute(
                data, Integral{.comp = c_x}, Min{.comp = c_x}, l3_norm,
                Max{.comp = c_x}, L2Norm{.comp = c_vector_0, .num_comps = 2});
            CHECK(together[0] ==
                  reducer.compute(data, Integral{.comp = c_x})[0]);
            CHECK(together[1] == reducer.compute(data, Min{.comp = c_x})[0]);
            CHECK(together[2] == reducer.compute(data, l3_norm)[0]);
            CHECK(together[3] == reducer.compute(data, Max{.comp = c_x})[0]);
            CHECK(together[4] ==
                  reducer.compute(
                      data, L2Norm{.comp = c_vector_0, .num_comps = 2})[0]);
        }
    }
    amrex::Finalize();
}
