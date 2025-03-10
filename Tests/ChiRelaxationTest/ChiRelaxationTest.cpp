/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */
// Doctest header
#include "doctest.h"

// Test header
#include "ChiRelaxationTest.hpp"

// Common test headers
#include "InitialData.hpp"
#include "doctestCLIArgs.hpp"

// GRTeclyn headers
#include "ChiRelaxation.hpp"
#include "DefaultPotential.hpp"
#include "ScalarField.hpp"

// AMReX headers
#include "AMReX.H"
#include "AMReX_FArrayBox.H"
#include "AMReX_MultiFab.H"

#ifdef AMREX_USE_HDF5
#include <AMReX_PlotFileUtilHDF5.H>
#endif

// System headers
#include <array>
#include <cstdlib>
#include <iostream>
#include <string>

void run_chi_relaxation_test()
{
    int amrex_argc    = doctest::cli_args.argc();
    char **amrex_argv = doctest::cli_args.argv();

    // NOLINTNEXTLINE(bugprone-casting-through-void)
    amrex::Initialize(amrex_argc, amrex_argv, true, MPI_COMM_WORLD);
    {

        constexpr int num_cells  = 32;
        constexpr int num_ghosts = 3;
        constexpr double dx      = 0.25 / (num_cells - 1);

        amrex::Box box(
            amrex::IntVect(0, 0, 0),
            amrex::IntVect(num_cells - 1, num_cells - 1, num_cells - 1));

        amrex::Box ghosted_box = box;
        ghosted_box.grow(num_ghosts);

        amrex::BoxArray box_array{box};
        amrex::RealVect dx_Vect{dx};
        amrex::RealBox real_box{box, dx_Vect.dataPtr(),
                                amrex::RealVect::Zero.dataPtr()};
        int coord_sys = 0;
        amrex::Geometry geom{box, &real_box, coord_sys};
        amrex::DistributionMapping distribution_mapping{box_array};
        amrex::MFInfo mf_info;
        mf_info.SetArena(amrex::The_Managed_Arena());

        amrex::MultiFab in_mf{box_array, distribution_mapping, NUM_VARS,
                              num_ghosts, mf_info};

        const auto &in_array = in_mf.arrays();

        // NOLINTBEGIN(bugprone-easily-swappable-parameters)
        amrex::ParallelFor(
            in_mf, in_mf.nGrowVect(),
            [=] AMREX_GPU_DEVICE(int ibox, int i, int j, int k)
            // NOLINTEND(bugprone-easily-swappable-parameters)
            {
                const amrex::IntVect iv{i, j, k};
                const amrex::RealVect coords = amrex::RealVect{iv} * dx;
                amrex::Real x                = coords[0];
                amrex::Real y                = coords[1];
                amrex::Real z                = coords[2];

                random_ccz4_initial_data(iv, in_array[ibox], coords);

                random_matter_bssn_initial_data(iv, in_array[ibox], coords);
            });

        using DefaultScalarField = ScalarField<DefaultPotential>;

        amrex::Real G_Newton    = 1.0;
        amrex::Real relax_speed = 0.1;

        ChiRelaxation<DefaultScalarField> chi_relaxation{
            DefaultScalarField(), dx, relax_speed, G_Newton};

        amrex::MultiFab out_mf{box_array, distribution_mapping, NUM_VARS, 0,
                               mf_info};

        const auto &in_c_array   = in_mf.const_arrays();
        const auto &out_mf_array = out_mf.arrays();

        amrex::ParallelFor(out_mf,
                           [=] AMREX_GPU_DEVICE(int ibox, int i, int j, int k)
                           {
                               chi_relaxation.compute(i, j, k,
                                                      out_mf_array[ibox],
                                                      in_c_array[ibox]);
                           });

        // GPU barrier
        amrex::Gpu::streamSynchronize();

        CHECK(!out_mf.contains_nan());

#if AMREX_USE_HDF5
        std::string grteclyn_hdf5_file = "ChiRelaxationTest/ChiRelaxationTest";

        // open the hdf5 file for writing
        amrex::WriteSingleLevelPlotfileHDF5(
            grteclyn_hdf5_file, out_mf, StateVariables::names, geom, 0.0, 0);

        std::cout.flush();

        std::string h5diff_tol         = "1e-10";
        std::string grchombo_hdf5_file = "ChiRelaxationTest/"
                                         "GRChomboChiRelaxationTest.h5";
        // the GRChombo comparison file is created in
        // GRChombo/Tests/ChiRelaxationTest in the branch
        // tests/MatterCCZ4RHS

        std::string hdf5_internal_path = "/level_0/data:datatype=0";

        std::string h5diff_command  = "h5diff";
        h5diff_command             += " -d " + h5diff_tol;
        h5diff_command +=
            " " + grteclyn_hdf5_file + ".h5" + " " + grchombo_hdf5_file;
        h5diff_command += " " + hdf5_internal_path + " " + hdf5_internal_path;

        INFO("Test command: " << h5diff_command);

        int h5diff_status = std::system(h5diff_command.c_str());
        int h5diff_retval = -1;

        // Use POSIX macros to get the exit code
        if (WIFEXITED(h5diff_status))
        {
            h5diff_retval = WEXITSTATUS(h5diff_status);
        }

        CHECK(h5diff_retval == 0);

#endif
    }
    amrex::Finalize();
}
