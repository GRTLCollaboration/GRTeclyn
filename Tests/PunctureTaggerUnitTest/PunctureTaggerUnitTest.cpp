/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// Doctest header
#include "doctest.h"

// Test header
#include "PunctureTaggerUnitTest.hpp"

// Common includes
#include "doctestCLIArgs.hpp"

// GRTeclyn includes
#include "GRParmParse.hpp"
#include "PunctureTagger.hpp"

// AMReX includes
#include <AMReX.H>
#include <AMReX_BaseFab.H>

#include <array>

void run_puncture_tagger_unit_test()
{
    int amrex_argc    = doctest::cli_args.argc();
    char **amrex_argv = doctest::cli_args.argv();
    // NOLINTNEXTLINE(bugprone-casting-through-void) // Open MPI triggers this
    amrex::Initialize(amrex_argc, amrex_argv);
    {
        GRParmParse pp;
        pp.add("puncture_tagging.level_separation", 1.5);
        pp.add("puncture_tagging.finest_level_factor", 2.0);

        constexpr int max_level = 6;
        using tagger_t          = PunctureTagger<2>;
        const std::array<amrex::Real, AMREX_SPACEDIM * 2> puncture_coords{};

        const auto check_max_levels =
            [&](const std::array<amrex::Real, 2> &a_masses,
                const std::array<int, 2> &a_expected_levels)
        {
            const tagger_t tagger(1.0, 0, max_level, puncture_coords, a_masses);
            CHECK(tagger.get_puncture_max_level(0) == a_expected_levels[0]);
            CHECK(tagger.get_puncture_max_level(1) == a_expected_levels[1]);
        };

        check_max_levels({1.0, 1.0}, {6, 6});
        check_max_levels({1.0, 1.999}, {6, 6});
        check_max_levels({1.0, 2.0}, {6, 5});
        check_max_levels({1.0, 3.999}, {6, 5});
        check_max_levels({1.0, 4.0}, {6, 4});
        check_max_levels({1.0, 8.0}, {6, 3});
        check_max_levels({1.0, 128.0}, {6, 0});
        check_max_levels({4.0, 1.0}, {4, 6});

        const std::array<amrex::Real, AMREX_SPACEDIM * 2>
            separated_puncture_coords{0.5, 0.5, 0.5, 16.5, 0.5, 0.5};
        const std::array<amrex::Real, 2> unequal_masses{1.0, 4.0};

        const auto is_tagged = [&](const int a_level, const int a_ix)
        {
            const amrex::IntVect cell(a_ix, 0, 0);
            amrex::BaseFab<amrex::TagBox::TagType> tags(amrex::Box(cell, cell),
                                                        1);
            tags.setVal(amrex::TagBox::CLEAR);
            const tagger_t tagger(1.0, a_level, max_level,
                                  separated_puncture_coords, unequal_masses);
            tagger(a_ix, 0, 0, tags.array());
            return tags(cell, 0) == amrex::TagBox::SET;
        };

        CHECK(is_tagged(5, 0));  // the smaller puncture reaches level 6
        CHECK(is_tagged(3, 16)); // the larger puncture reaches level 4
        CHECK_FALSE(is_tagged(4, 16));
    }
    amrex::Finalize();
}
