/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// AMReX includes
#include <AMReX.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_FabArray.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>

#ifdef AMREX_USE_HDF5
#include <AMReX_PlotFileUtilHDF5.H>
#endif

// Doctest includes
#include "doctest.h"
#include "doctestCLIArgs.hpp"

#include <iomanip>
#include <iostream>
#include <sys/time.h>

// common includes
#include "InitialData.hpp" //includes StateVariables.hpp

// test header
#include "EMTensorTest.hpp"

// GRTeclyn includes
#include "DefaultPotential.hpp"
#include "EMTensor.hpp"
#include "ScalarField.hpp"
#include "simd.hpp"
#include <array>

void run_emtensor_test()
{
    int amrex_argc    = doctest::cli_args.argc();
    char **amrex_argv = doctest::cli_args.argv();

    // NOLINTNEXTLINE(bugprone-casting-through-void) // Open MPI triggers this
    amrex::Initialize(amrex_argc, amrex_argv, true, MPI_COMM_WORLD);
    {

        constexpr int num_cells  = 32;
        constexpr int num_ghosts = 3;
        constexpr amrex::Real dx = 0.25 / (num_cells - 1);
        amrex::Box box{amrex::IntVect::TheZeroVector(),
                       amrex::IntVect{num_cells - 1}};
        auto ghosted_box = box;
        ghosted_box.grow(num_ghosts);

        amrex::BoxArray box_array{box};

        amrex::RealVect dx_Vect{dx};
        amrex::RealBox real_box{box, dx_Vect.dataPtr(),
                                amrex::RealVect::Zero.dataPtr()};

        int coord_sys = 0; // Cartesian

        amrex::Geometry geom{box, &real_box, coord_sys};
        amrex::DistributionMapping distribution_mapping{box_array};
        amrex::MFInfo mf_info;
        mf_info.SetArena(amrex::The_Managed_Arena());

        amrex::MultiFab in_mf{box_array, distribution_mapping, NUM_VARS,
                              num_ghosts, mf_info};

        const auto &in_arrays = in_mf.arrays();

        // NOLINTBEGIN(bugprone-easily-swappable-parameters)
        amrex::ParallelFor(
            in_mf, in_mf.nGrowVect(),
            [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
            // NOLINTEND(bugprone-easily-swappable-parameters)
            {
                const amrex::IntVect iv{i, j, k};
                const amrex::RealVect coords = amrex::RealVect{iv} * dx;
                amrex::Real x                = coords[0];
                amrex::Real y                = coords[1];
                amrex::Real z                = coords[2];

                random_ccz4_initial_data(iv, in_arrays[box_no], coords);

                random_matter_bssn_initial_data(iv, in_arrays[box_no], coords);
            });

        amrex::Gpu::streamSynchronize();

        // Setup scalar field calculations
        using DefaultScalarField = ScalarField<DefaultPotential>;

        // // Constructor for EMTensor
        constexpr int dcomp_rho = 0;

        constexpr int num_comps = 1;

        amrex::MultiFab out_mf{box_array, distribution_mapping, num_comps,
                               num_ghosts, mf_info};

        const auto &in_c_arrays = in_mf.const_arrays();

        double time = 0.0;
        int *bcrec  = nullptr;
        int level   = 0;

        EMTensor<DefaultScalarField,
                 EMTensorOptions::justEnergyDensity>::compute_mf(out_mf,
                                                                 dcomp_rho,
                                                                 num_comps,
                                                                 in_mf, geom,
                                                                 time, bcrec,
                                                                 level);

#if AMREX_USE_HDF5
        amrex::Vector<std::string> var_names = {"rho"};

        std::string grteclyn_hdf5_file = "EMTensorTest/EMTensorTestOut";

        // open the hdf5 file for writing
        amrex::WriteSingleLevelPlotfileHDF5(grteclyn_hdf5_file, out_mf,
                                            var_names, geom, 0.0, 0);

        std::cout.flush();

        std::string h5diff_tol         = "1e-10";
        std::string grchombo_hdf5_file = "EMTensorTest/GRChomboEMTensorTest.h5";

        std::string hdf5_internal_path = "/level_0/data:datatype=0";

        std::string hdf5_tool      = "h5diff";
        std::string h5diff_command = hdf5_tool + " -d " + h5diff_tol;

        h5diff_command = h5diff_command + " " + grteclyn_hdf5_file + ".h5 " +
                         grchombo_hdf5_file + " " + hdf5_internal_path;

        INFO("Test command : " << h5diff_command);

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
