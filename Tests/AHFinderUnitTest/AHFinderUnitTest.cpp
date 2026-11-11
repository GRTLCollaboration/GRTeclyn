// Test include
#include "AHFinderUnitTest.hpp"

// Doctest header
#include "doctest.h"

// Common includes
#include "doctestCLIArgs.hpp"

// AMReX includes
#include "AMReX.H"
#include <AMReX_Print.H>

// Base includes
#include "BHAmr.hpp"
#include "DefaultLevelBld.hpp"
#include "GRParmParse.hpp"

// Problem specific includes
#include "AHFinder.hpp"
#include "AHFinderLevel.hpp"
#include "SimulationParameters.hpp"

// Others
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <vector>

void run_ah_finder_unit_test()
{
    // Use an input file that is in the same directory as this file for the
    // second argument
    std::filesystem::path this_file(__FILE__);
    std::filesystem::path input_file =
        this_file.parent_path() / std::filesystem::path("params_test.txt");
    char *input_file_c_str = strdup(input_file.c_str());

    auto new_args = doctest::cli_args;
    new_args.insert(1, input_file_c_str);

    int new_argc    = new_args.argc();
    char **new_argv = new_args.argv();

    // NOLINTNEXTLINE(bugprone-casting-through-void) // Open MPI triggers this
    amrex::Initialize(
        new_argc, new_argv,
        std::function<void()>(SimulationParameters::check_params));
    {
        // Simulation parameters
        GRParmParse pp;
        AHFinderLevel::variableSetUp();

        // Set up the AMR object. BHAmr owns the AH finders.
        constexpr int num_punctures = 2;
        DefaultLevelBld<AHFinderLevel> ahfinder_test_level_fact;
        BHAmr<num_punctures> bh_amr(&ahfinder_test_level_fact);

        amrex::Real stop_time{};
        pp.get("evolution.stop_time", stop_time);
        bh_amr.init(0., stop_time);

        // BHAmr creates one finder per puncture from the ah_finder
        // parameters. Find each individual horizon on the initial data, then
        // again after each of a few evolution steps, each find() starting
        // from the last converged surface. The gauge moves the horizons in
        // coordinates, but their areas, and so their irreducible masses,
        // should barely change over this short evolution.
        REQUIRE(bh_amr.ah_finder_enabled());

        std::array<std::vector<double>, num_punctures> h_initial;
        std::array<amrex::Real, num_punctures> mass_initial{};

        for (int ibh = 0; ibh < num_punctures; ++ibh)
        {
            AHFinder &finder = bh_amr.get_ah_finder(ibh);
            CHECK(finder.find(bh_amr.cumTime(), bh_amr.dtLevel(0)));
            h_initial[ibh]    = finder.h();
            mass_initial[ibh] = finder.mass();
        }

        // The data is symmetric under reflection in x, so both horizons
        // should have the same mass
        CHECK(mass_initial[1] == doctest::Approx(mass_initial[0]));

        int max_steps{};
        pp.get("evolution.max_steps", max_steps);

        // Relative tolerance on the change in the irreducible mass
        constexpr double mass_tolerance = 1e-3;

        for (int step = 0; step < max_steps; ++step)
        {
            bh_amr.coarseTimeStep(stop_time);

            for (int ibh = 0; ibh < num_punctures; ++ibh)
            {
                AHFinder &finder = bh_amr.get_ah_finder(ibh);
                CHECK(finder.find(bh_amr.cumTime(), bh_amr.dtLevel(0)));
                CHECK(
                    finder.mass() ==
                    doctest::Approx(mass_initial[ibh]).epsilon(mass_tolerance));
            }
        }

        for (int ibh = 0; ibh < num_punctures; ++ibh)
        {
            const std::vector<double> &h = bh_amr.get_ah_finder(ibh).h();

            double max_h_change = 0.0;
            for (std::size_t i = 0; i < h.size(); ++i)
            {
                max_h_change =
                    std::max(max_h_change, std::abs(h[i] - h_initial[ibh][i]));
            }
            amrex::Print() << "Max change in h of horizon " << ibh + 1
                           << " over " << max_steps
                           << " steps = " << max_h_change << "\n";
            // The data has genuinely evolved, so the surface should have moved
            CHECK(max_h_change > 0.0);
        }
    }

    amrex::Finalize();
}
