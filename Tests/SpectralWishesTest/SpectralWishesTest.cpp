/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// Doctest header
#include "doctest.h"

// Test header
#include "SpectralWishesTest.hpp"
#include "SpectralWishesTestLevel.hpp"

// Common test headers
// #include "InitialData.hpp"
#include "doctestCLIArgs.hpp"

// GRTeclyn headers
#include "GRAmr.hpp"

// AMReX headers
#include "AMReX.H"
#include "AMReX_FArrayBox.H"
#include "AMReX_MultiFab.H"

// System headers
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>

void run_spectral_wishes_test()
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

    // int amrex_argc    = doctest::cli_args.argc();
    // char **amrex_argv = doctest::cli_args.argv();

    // NOLINTNEXTLINE(bugprone-casting-through-void) // Open MPI triggers this
    amrex::Initialize(
        new_argc, new_argv,
        std::function<void()>(SimulationParameters::check_params));
    {
        GRParmParse pp; // NOLINT(readability-identifier-length)

        DefaultLevelBld<SpectralWishesTestLevel> level_factory;

        GRAmr gr_amr(&level_factory);
        double stop_time{};
        pp.get("evolution.stop_time", stop_time);
        gr_amr.init(0., stop_time);

        int max_steps{};
        pp.get("evolution.max_steps", max_steps);

        while ((gr_amr.okToContinue() != 0) &&
               (gr_amr.levelSteps(0) < max_steps || max_steps < 0) &&
               (gr_amr.cumTime() < stop_time || stop_time < 0.0))
        {
            gr_amr.coarseTimeStep(stop_time);
        }
    }
    amrex::Finalize();
    // constexpr int num_cells  = 32;
    // constexpr int num_ghosts = 3;
    // constexpr amrex::Real dx = 0.5 / (num_cells - 1);

    // amrex::Box box(
    //     amrex::IntVect(0, 0, 0),
    //     amrex::IntVect(num_cells - 1, num_cells - 1, num_cells - 1));

    // amrex::Box ghosted_box = box;
    // ghosted_box.grow(num_ghosts);

    // amrex::BoxArray box_array{box};
    // amrex::RealVect dx_Vect{dx};
    // amrex::RealBox real_box{box, dx_Vect.dataPtr(),
    //                         amrex::RealVect::Zero.dataPtr()};
    // int coord_sys = 0;
    // amrex::Geometry geom{box, &real_box, coord_sys};
    // amrex::DistributionMapping distribution_mapping{box_array};
    // amrex::MFInfo mf_info;
    // mf_info.SetArena(amrex::The_Managed_Arena());

    // amrex::MultiFab in_mf{box_array, distribution_mapping, NUM_VARS,
    //                       num_ghosts, mf_info};
    // in_mf.setVal(0.0); // initialise to zero

    // const auto &in_array = in_mf.arrays();

    // // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    // amrex::ParallelFor(
    //     in_mf, in_mf.nGrowVect(),
    //     [=] AMREX_GPU_DEVICE(int ibox, int i, int j, int k)
    //     // NOLINTEND(bugprone-easily-swappable-parameters)
    //     {
    //         const amrex::IntVect iv{i, j, k};
    //         const amrex::RealVect coords = amrex::RealVect{iv} * dx;
    //         amrex::Real x                = coords[0];
    //         amrex::Real y                = coords[1];
    //         amrex::Real z                = coords[2];

    //         random_ccz4_initial_data(iv, in_array[ibox], coords);

    //         random_matter_bssn_initial_data(iv, in_array[ibox], coords);
    //     });

    // constexpr int num_bssn_matter_vars = c_Pi + 1;
    // constexpr int dcomp                = 0;

    // int num_comp_constraints = 1 + AMREX_SPACEDIM; // ham + moms

    // // out_mf only contains constraints
    // amrex::MultiFab out_mf{box_array, distribution_mapping,
    //                        num_comp_constraints, 0, mf_info};

    // amrex::FArrayBox out_fab{box, num_comp_constraints,
    //                          amrex::The_Managed_Arena()};

    // // const auto &in_c_array    = in_mf.const_arrays();
    // const auto &out_mf_array  = out_mf.arrays();
    // const auto &out_fab_array = out_fab.array();

    // CHECK(!out_mf.contains_nan());
}
