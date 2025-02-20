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
#include "BHAMR.hpp"
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
    int amrex_argc    = doctest::cli_args.argc();
    char **amrex_argv = doctest::cli_args.argv();

    // Add an inputs file to our arguments
    int new_argc = amrex_argc + 1;
    std::vector<char *> new_args(new_argc);
    new_args[0] = amrex_argv[0];

    // Use an input file that is in the same directory as this file for the
    // second argument
    std::filesystem::path this_file(__FILE__);
    std::filesystem::path input_file =
        this_file.parent_path() / std::filesystem::path("test.inputs");
    char input_file_c_str[PATH_MAX];
    std::strcpy(input_file_c_str, input_file.c_str());
    new_args[1] = input_file_c_str;

    for (int iarg = 2; iarg < new_argc; ++iarg)
    {
        new_args[iarg] = amrex_argv[iarg - 1];
    }
    char **new_argv = new_args.data();

    // NOLINTNEXTLINE(bugprone-casting-through-void) // Open MPI triggers this
    amrex::Initialize(new_argc, new_argv, true, MPI_COMM_WORLD);
    {
        GRParmParse pp; // NOLINT(readability-identifier-length)
        SimulationParameters sim_params(pp);

        GRAMR::set_simulation_parameters(sim_params);

        DefaultLevelFactory<PunctureTrackerLevel> level_factory;

        BHAMR<2> bh_amr(&level_factory);
        bh_amr.init(0., sim_params.stop_time);

        while ((bh_amr.okToContinue() != 0) &&
               (bh_amr.levelSteps(0) < sim_params.max_steps ||
                sim_params.max_steps < 0) &&
               (bh_amr.cumTime() < sim_params.stop_time ||
                sim_params.stop_time < 0.0))
        {
            bh_amr.coarseTimeStep(sim_params.stop_time);
        }
    }
    amrex::Finalize();
}
