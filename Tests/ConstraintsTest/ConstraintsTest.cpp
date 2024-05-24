/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// Doctest header
#include "doctest.h"

// Test header
#include "ConstraintsTest.hpp"

// Common Test headers
#include "InitialData.hpp"
#include "doctestCLIArgs.hpp"

// Our includes
#include "NewConstraints.hpp"

// AMReX headers
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#ifdef AMREX_USE_HDF5
#include <AMReX_PlotFileUtilHDF5.H>
#endif

// System headers
#include <array>
#include <cstdlib>
#include <iostream>
#include <string>

void run_constraints_test()
{
    int amrex_argc    = doctest::cli_args.argc();
    char **amrex_argv = doctest::cli_args.argv();
    amrex::Initialize(amrex_argc, amrex_argv, true, MPI_COMM_WORLD);
    {
        // Set up grid
        constexpr int num_cells  = 16;
        constexpr int num_ghosts = 3;
        constexpr amrex::Real dx = 1.0 / (num_cells - 1);

        amrex::Box box{amrex::IntVect::TheZeroVector(),
                       amrex::IntVect{num_cells - 1}};
        auto ghosted_box = box;
        ghosted_box.grow(num_ghosts);

        amrex::BoxArray box_array{box};
        amrex::DistributionMapping distribution_mapping{box_array};
        amrex::MFInfo mf_info;
        mf_info.SetArena(amrex::The_Managed_Arena());

        amrex::MultiFab in_mf{box_array, distribution_mapping, NUM_CCZ4_VARS,
                              num_ghosts, mf_info};

        // Calculate initial data
        const auto &in_arrays = in_mf.arrays();
        amrex::ParallelFor(
            in_mf, in_mf.nGrowVect(),
            [=] AMREX_GPU_DEVICE(int ibox, int i, int j, int k)
            {
                const amrex::IntVect iv{i, j, k};
                const amrex::RealVect coords = amrex::RealVect{iv} * dx;

                random_ccz4_initial_data(iv, in_arrays[ibox], coords);
            });

        amrex::Gpu::streamSynchronize();

        constexpr int dcomp = 0;
        int iham            = dcomp;
        Interval imom       = Interval(dcomp + 1, dcomp + AMREX_SPACEDIM);
        Constraints constraints(dx, iham, imom);

        int num_constraints_comps = 4;
        int num_out_ghosts        = 0;
        amrex::MultiFab out_mf{box_array, distribution_mapping,
                               num_constraints_comps, num_out_ghosts, mf_info};
        const auto &in_c_arrays = in_mf.const_arrays();
        const auto &out_arrays  = out_mf.arrays();
        amrex::ParallelFor(out_mf,
                           [=] AMREX_GPU_DEVICE(int ibox, int i, int j, int k) {
                               constraints.compute(i, j, k, out_arrays[ibox],
                                                   in_c_arrays[ibox]);
                           });

        amrex::Gpu::streamSynchronize();

        // Write to HDF5 plot file
#ifdef AMREX_USE_HDF5
        amrex::RealVect dx_Vect{dx};
        amrex::RealBox real_box{box, dx_Vect.dataPtr(),
                                amrex::RealVect::Zero.dataPtr()};

        int coord_sys = 0; // Cartesian

        amrex::Geometry geom{box, &real_box, coord_sys};

        std::string this_test_dir = "ConstraintsTest/";
        std::string hdf5_out_stem = this_test_dir + "ConstraintsOut";

        amrex::WriteSingleLevelPlotfileHDF5(
            hdf5_out_stem, out_mf, Constraints::var_names, geom, 1.0, 0);

        // Apparently this is necessary before calling std::system if h5diff
        // writes to the screen
        std::cout.flush();

        std::string h5diff_toll        = "1.0e-10";
        std::string grteclyn_hdf5_file = hdf5_out_stem + ".h5";
        std::string grchombo_hdf5_file =
            this_test_dir + "ConstraintsGRChombo.hdf5";
        std::string hdf5_internal_path = "/level_0/data:datatype=0";

        // Let's hope this is in our PATH if we're building with HDF5
        std::string h5diff_command  = "h5diff";
        h5diff_command             += " -d " + h5diff_toll;
        h5diff_command += " " + grteclyn_hdf5_file + " " + grchombo_hdf5_file;
        h5diff_command += " " + hdf5_internal_path + " " + hdf5_internal_path;
        INFO("h5diff command: " << h5diff_command);
        INFO("Run command manually with -r flag to print difference");

        // In an ideal world, we wouldn't rely on running an external program
        // but I can't think of a simple way to do this nicely.
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