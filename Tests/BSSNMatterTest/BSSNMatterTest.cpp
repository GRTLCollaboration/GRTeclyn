/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */
// Doctest header
#include "doctest.h"

// Test header
#include "BSSNMatterTest.hpp"

// Common test headers
#include "InitialData.hpp"
#include "doctestCLIArgs.hpp"

// GRTeclyn headers
#include "DefaultPotential.hpp"
#include "FourthOrderDerivatives.hpp"
#include "MatterCCZ4RHS.hpp"
#include "NewMatterConstraints.hpp"
#include "ScalarField.hpp"

// AMReX headers
#include "AMReX.H"
#include "AMReX_FArrayBox.H"
#include "AMReX_MultiFab.H"

#ifdef AMREX_USE_HDF5
#include "H5Cpp.h"
#include <AMReX_PlotFileUtilHDF5.H>
using namespace H5;
#endif

// System headers
#include <array>
#include <cstdlib>
#include <iostream>
#include <string>

void run_bssn_matter_test()
{
    int amrex_argc    = doctest::cli_args.argc();
    char **amrex_argv = doctest::cli_args.argv();
    amrex::Initialize(amrex_argc, amrex_argv, true, MPI_COMM_WORLD);
    {

        constexpr int num_cells  = 32;
        constexpr int num_ghosts = 3;
        constexpr double dx      = 0.5 / (num_cells - 1);

        amrex::Box box(
            amrex::IntVect(0, 0, 0),
            amrex::IntVect(num_cells - 1, num_cells - 1, num_cells - 1));

        amrex::Box ghosted_box = box;
        ghosted_box.grow(num_ghosts);

        amrex::BoxArray box_array{box};
        amrex::DistributionMapping distribution_mapping{box_array};
        amrex::MFInfo mf_info;
        mf_info.SetArena(amrex::The_Managed_Arena());

        amrex::MultiFab in_mf{box_array, distribution_mapping, NUM_VARS,
                              num_ghosts, mf_info};

        amrex::FArrayBox in_fab{ghosted_box, NUM_VARS,
                                amrex::The_Managed_Arena()};

        const auto &in_array     = in_mf.arrays();
        const auto &in_fab_array = in_fab.array();

        amrex::ParallelFor(
            in_mf, in_mf.nGrowVect(),
            [=] AMREX_GPU_DEVICE(int ibox, int i, int j, int k)
            {
                const amrex::IntVect iv{i, j, k};
                const amrex::RealVect coords = amrex::RealVect{iv} * dx;
                amrex::Real x                = coords[0];
                amrex::Real y                = coords[1];
                amrex::Real z                = coords[2];

                random_ccz4_initial_data(iv, in_array[ibox], coords);

                // the initial data doesn't include phi or Pi so do it here:
                in_array[ibox](i, j, k, c_phi) =
                    0.34578 + 0.26898 * x + 0.54348 * x * x +
                    0.33487 * x * y * y * y + 0.79469 * y * z +
                    0.30515 * z * z + 1.88385 * z * z * z * z;
                in_array[ibox](i, j, k, c_Pi) =
                    0.65668 + 0.20188 * x + 0.34348 * x * x +
                    0.31787 * x * y * y * y + 0.88469 * y * z +
                    0.10515 * z * z + 1.88385 * z * z * z * z;
            });

        amrex::ParallelFor(
            ghosted_box,
            [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                const amrex::IntVect iv{i, j, k};
                const amrex::RealVect coords = amrex::RealVect{iv} * dx;
                amrex::Real x                = coords[0];
                amrex::Real y                = coords[1];
                amrex::Real z                = coords[2];

                random_ccz4_initial_data(iv, in_fab_array, coords);

                // Theta is set to zero in BSSN
                in_fab_array(i, j, k, c_Theta) = 0.0;

                // the initial data doesn't include phi or Pi so do it here:
                in_fab_array(i, j, k, c_phi) =
                    0.34578 + 0.26898 * x + 0.54348 * x * x +
                    0.33487 * x * y * y * y + 0.79469 * y * z +
                    0.30515 * z * z + 1.88385 * z * z * z * z;
                in_fab_array(i, j, k, c_Pi) =
                    0.65668 + 0.20188 * x + 0.34348 * x * x +
                    0.31787 * x * y * y * y + 0.88469 * y * z +
                    0.10515 * z * z + 1.88385 * z * z * z * z;
            });

        CCZ4_params_t<MovingPunctureGauge::params_t> ccz4_params;
        ccz4_params.kappa1            = 0.0;
        ccz4_params.kappa2            = 0.0;
        ccz4_params.kappa3            = 0.0;
        ccz4_params.shift_Gamma_coeff = 0.75;
        ccz4_params.lapse_advec_coeff = 1.0;
        ccz4_params.lapse_power       = 1.0;
        ccz4_params.lapse_coeff       = 2.0;
        ccz4_params.shift_advec_coeff = 0.0;
        ccz4_params.eta               = 1.0;

        amrex::Real sigma = 0.1;

        using DefaultScalarPotential = ScalarField<DefaultPotential>;

        double G_Newton = 1.0;

        MatterCCZ4RHS<DefaultScalarPotential, MovingPunctureGauge,
                      FourthOrderDerivatives>
            current_ccz4_rhs{DefaultScalarPotential(DefaultPotential()),
                             ccz4_params,
                             dx,
                             sigma,
                             CCZ4RHS<>::USE_BSSN,
                             G_Newton};

        // Set up the constraints
        constexpr int dcomp = NUM_VARS;
        int iham            = dcomp;
        Interval imom       = Interval(dcomp + 1, dcomp + AMREX_SPACEDIM);
        MatterConstraints<DefaultScalarPotential> constraints(
            DefaultScalarPotential(DefaultPotential()), dx, G_Newton, iham,
            imom);

        int num_constraints_comps = 4;

        amrex::MultiFab out_mf{box_array, distribution_mapping, NUM_VARS, 0,
                               mf_info};

        amrex::FArrayBox out_fab{box, NUM_VARS + num_constraints_comps,
                                 amrex::The_Managed_Arena()};

        const auto &in_c_array    = in_mf.const_arrays();
        const auto &out_mf_array  = out_mf.arrays();
        const auto &out_fab_array = out_fab.array();

        // Do the CCZ4RHS and constraints calculation in the same loop

        const auto &in_fab_c_array = in_fab.const_array();

        amrex::ParallelFor(box,
                           [=] AMREX_GPU_DEVICE(int i, int j, int k)
                           {
                               current_ccz4_rhs.compute(i, j, k, out_fab_array,
                                                        in_fab_c_array);

                               constraints.compute(i, j, k, out_fab_array,
                                                   in_fab_c_array);
                           });

        amrex::ParallelFor(out_mf,
                           [=] AMREX_GPU_DEVICE(int ibox, int i, int j, int k) {
                               current_ccz4_rhs.compute(i, j, k,
                                                        out_mf_array[ibox],
                                                        in_c_array[ibox]);
                           });

        // GPU barrier
        amrex::Gpu::streamSynchronize();

        if (out_mf.contains_nan())
        {
            amrex::Abort("NaN found in output multifab");
            CHECK(false);
        }
        else
        {
            CHECK(true);
        }

#if AMREX_USE_HDF5
        int coord_sys = 0;
        amrex::Vector<std::string> var_names;
        for (int i = 0; i < NUM_VARS; i++)
            var_names.push_back(StateVariables::names[i]);
        for (int i = 0; i < num_constraints_comps; i++)
            var_names.push_back(Constraints::var_names[i]);

        const H5std_string grteclyn_hdf5_file =
            "BSSNMatterTest/BSSNMatterTest.h5";

        // open the hdf5 file for writing
        hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
#if AMREX_USE_MPI
        MPI_Info mpi_info = MPI_INFO_NULL;
        H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, mpi_info);
#endif
        hid_t fid =
            H5Fcreate(grteclyn_hdf5_file.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                      plist_id); // H5F_ACC_TRUNC = if file exists open with
                                 // read/write access, otherwise create file

        // // create the group
        char level_name[8] = "level_0"; // only the one level

        hid_t grp =
            H5Gcreate(fid, level_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        for (int i = 0; i < NUM_VARS + num_constraints_comps; i++)
        {

            // // create the dataset
            char dataname[32] = "data:datatype=0";
            hsize_t hs_procsize[1]; // only 1 dimension
            hs_procsize[0] = box.length(0) * box.length(1) * box.length(2);
            //                         NUM_VARS; // total number of values
            //                         printed
            hid_t dataspace    = H5Screate_simple(1, hs_procsize, NULL);
            hid_t memdataspace = H5Screate_simple(1, hs_procsize, NULL);

            hid_t dcpl_id  = H5Pcreate(H5P_DATASET_CREATE);
            hid_t dxpl_col = H5Pcreate(H5P_DATASET_XFER);

            hid_t dataset =
                H5Dcreate(grp, var_names[i].c_str(), H5T_NATIVE_DOUBLE,
                          dataspace, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
            hid_t ret = H5Dwrite(
                dataset, H5T_NATIVE_DOUBLE, memdataspace, dataspace, dxpl_col,
                out_fab.dataPtr(i)); // could also pass in a parameter to
                                     // dataPtr e.g. dataPtr(c_Theta)
            H5Dclose(dataset);
            H5Sclose(memdataspace);
            H5Sclose(dataspace);
            H5Pclose(dxpl_col);
        }
        H5Gclose(grp);
        H5Pclose(plist_id);
        H5Fclose(fid);

        std::cout.flush();

        std::string h5diff_tol         = "1e-10";
        std::string grchombo_hdf5_file = "BSSNMatterTest/"
                                         "GRChomboBSSNMatterTest.h5";
        // the GRChombo comparison file is created in
        // GRChombo/Tests/MatterCCZ4Test

        // could also give an internal path
        // but leaving it off does the full file
        std::string hdf5_internal_path = "/level_0/";

        std::string h5diff_command  = "h5diff";
        h5diff_command             += " -d " + h5diff_tol;
        h5diff_command += " " + grteclyn_hdf5_file + " " + grchombo_hdf5_file;
        h5diff_command += " " + hdf5_internal_path + " " + hdf5_internal_path;

        amrex::Print() << "Test command: " << h5diff_command << std::endl;

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
