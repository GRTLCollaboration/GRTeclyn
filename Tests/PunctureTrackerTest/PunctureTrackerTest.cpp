/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// Doctest header
#include "doctest.h"

// Common test includes
#include "SimulationParameters.hpp"
#include "doctestCLIArgs.hpp"

// GRTeclyn includes
#include "BHAmr.hpp"
#include "DefaultLevelFactory.hpp"
#include "GRParmParse.hpp"
#include "SetupFunctions.hpp"

// Problem specific includes:
#include "PunctureTrackerLevel.hpp"

// System includes
#include <chrono>
#include <climits> // for PATH_MAX
#include <filesystem>
#include <iostream>

// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
void run_puncture_tracker_test()
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
        GRParmParse pp; // NOLINT(readability-identifier-length)

        DefaultLevelFactory<PunctureTrackerLevel> level_factory;

        BHAmr<2> bh_amr(&level_factory);
        double stop_time{};
        pp.get("evolution.stop_time", stop_time);
        bh_amr.init(0., stop_time);

        int max_steps{};
        pp.get("evolution.max_steps", max_steps);

        while ((bh_amr.okToContinue() != 0) &&
               (bh_amr.levelSteps(0) < max_steps || max_steps < 0) &&
               (bh_amr.cumTime() < stop_time || stop_time < 0.0))
        {
            bh_amr.coarseTimeStep(stop_time);
        }
    }
    amrex::Finalize();
}
