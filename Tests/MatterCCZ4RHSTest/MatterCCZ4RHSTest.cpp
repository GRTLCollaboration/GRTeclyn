/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */
// Doctest header
#include "doctest.h"

// Test header
#include "MatterCCZ4RHSTest.hpp"

// Common test headers
#include "InitialData.hpp"
#include "StateVariables.hpp"
// #include "CCZ4UserVariables.hpp"
#include "doctestCLIArgs.hpp"

// GRTeclyn headers
#include "DefaultPotential.hpp"
#include "FourthOrderDerivatives.hpp"
#include "MatterCCZ4RHS.hpp"
#include "NewMatterConstraints.hpp"
#include "ScalarField.hpp"

// Old GRTeclyn headers for comparison
#include "CCZ4RHS-fdf5a7a.hpp"
#include "FourthOrderDerivatives-fdf5a7a.hpp"

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

void run_matter_ccz4_rhs_test()
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
            in_mf,
            [=] AMREX_GPU_DEVICE(int ibox, int i, int j, int k)
            {
                const amrex::IntVect iv{i, j, k};
                const amrex::RealVect coords = amrex::RealVect{iv} * dx;

                random_ccz4_initial_data(iv, in_array[ibox], coords);
            });

        amrex::ParallelFor(
            ghosted_box,
            [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                const amrex::IntVect iv{i, j, k};
                const amrex::RealVect coords = amrex::RealVect{iv} * dx;

                random_ccz4_initial_data(iv, in_fab_array, coords);
            });

        CCZ4_params_t<MovingPunctureGauge::params_t> current_ccz4_params;
        current_ccz4_params.kappa1            = 0.1;
        current_ccz4_params.kappa2            = 0;
        current_ccz4_params.kappa3            = 1;
        current_ccz4_params.covariantZ4       = true;
        current_ccz4_params.lapse_advec_coeff = 0.0;
        current_ccz4_params.lapse_power       = 1.0;
        current_ccz4_params.lapse_coeff       = 2.0;
        current_ccz4_params.shift_Gamma_coeff = 0.75;
        current_ccz4_params.shift_advec_coeff = 0;
        current_ccz4_params.eta               = 1.82;

        amrex::Real sigma = 0.3;

        using DefaultScalarPotential = ScalarField<DefaultPotential>;

        double G_Newton = 1.0;

        MatterCCZ4RHS<DefaultScalarPotential, MovingPunctureGauge,
                      FourthOrderDerivatives>
            current_ccz4_rhs{DefaultScalarPotential(DefaultPotential()),
                             current_ccz4_params,
                             dx,
                             sigma,
                             CCZ4RHS<>::USE_CCZ4,
                             G_Newton};

        amrex::MultiFab out_mf{box_array, distribution_mapping, NUM_VARS, 0,
                               mf_info};

        amrex::FArrayBox out_fab{box, NUM_VARS, amrex::The_Managed_Arena()};
        amrex::FArrayBox diff_fab{box, NUM_VARS, amrex::The_Managed_Arena()};

        const auto &in_c_array    = in_mf.const_arrays();
        const auto &out_mf_array  = out_mf.arrays();
        const auto &out_fab_array = out_fab.array();
        const auto &diff_array    = diff_fab.array();

        // Do the current and old CCZ4RHS calculation in the same loop

        const auto &in_fab_c_array = in_fab.const_array();

        amrex::ParallelFor(box,
                           [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                               current_ccz4_rhs.compute(i, j, k, out_fab_array,
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
            CHECK(true);

        if (out_fab.contains_nan())
        {
            amrex::Abort("NaN found in output FArray");
            CHECK(false);
        }
        else
            CHECK(true);

#if AMREX_USE_HDF5
        int coord_sys                        = 0;
        amrex::Vector<std::string> var_names = {"Theta"};
        amrex::RealVect dx_Vect{dx};
        amrex::RealBox real_box{box, dx_Vect.dataPtr(),
                                amrex::RealVect::Zero.dataPtr()};
        amrex::Geometry geom{box, &real_box, coord_sys};
        // amrex::WriteSingleLevelPlotfileHDF5(
        //     "MatterCCZ4RHSTest/MatterCCZ4RHSOut", out_mf, var_names,
        //     geom, 0.0, 0);

        const H5std_string filename = "MatterCCZ4RHSTest/MatterCCZ4Out.h5";

        // open the hdf5 file for writing
        hid_t fid =
            H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                      H5P_DEFAULT); // H5F_ACC_TRUNC = if file exists open with
                                    // read/write access, otherwise create file

        // // create the group
        char level_name[8] = "level_0"; // only the one level
        hid_t grp =
            H5Gcreate(fid, level_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        // // create the dataset
        char dataname[64] = "data:datatype=0";
        hsize_t hs_procsize[1]; // only 1 dimension
        hs_procsize[0] = box.length(0) * box.length(1) * box.length(2) *
                         NUM_VARS; // total number of values printed
        hid_t dataspace    = H5Screate_simple(1, hs_procsize, NULL);
        hid_t memdataspace = H5Screate_simple(1, hs_procsize, NULL);

        hid_t dcpl_id  = H5Pcreate(H5P_DATASET_CREATE);
        hid_t dxpl_col = H5Pcreate(H5P_DATASET_XFER);

        hid_t dataset = H5Dcreate(grp, dataname, H5T_NATIVE_DOUBLE, dataspace,
                                  H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
        hid_t ret     = H5Dwrite(
            dataset, H5T_NATIVE_DOUBLE, memdataspace, dataspace, dxpl_col,
            out_fab.dataPtr()); // could also pass in a parameter to dataPtr
                                    // e.g. dataPtr(c_Theta)
        H5Dclose(dataset);
        H5Sclose(memdataspace);
        H5Sclose(dataspace);
        H5Pclose(dxpl_col);
        H5Gclose(grp);
        H5Fclose(fid);

        std::cout.flush();

        std::string h5diff_tol         = "1e-10";
        std::string grteclyn_hdf5_file = "MatterCCZ4RHSTest/MatterCCZ4Out.h5";
        std::string grchombo_hdf5_file =
            "MatterCCZ4RHSTest/GRChomboMatterCCZ4Out.h5"; // this file is
                                                          // created in
                                                          // GRChombo/Tests/MatterCCZ4Test
        std::string hdf5_internal_path = "/level_0/data:datatype=0";

        std::string h5diff_command  = "h5diff";
        h5diff_command             += " -d " + h5diff_tol;
        h5diff_command += " " + grteclyn_hdf5_file + " " + grchombo_hdf5_file;
        h5diff_command += " " + hdf5_internal_path + " " + hdf5_internal_path;

        int h5diff_retval = std::system(h5diff_command.c_str());

        CHECK(h5diff_retval == 0);

#endif

        amrex::Real max_diff = 0.0;
        amrex::IntVect max_diff_index{};

        // const int cout_precision = 17;
        // for (int ivar = 0; ivar < NUM_VARS; ++ivar)
        // {
        //     diff_fab.maxIndex<amrex::RunOn::Device>(box, max_diff,
        //                                             max_diff_index, ivar);

        //     INFO("Max diff for var "
        //          << UserVariables::variable_names[ivar] << ": "
        //          << std::setprecision(cout_precision) << max_diff << " at "
        //          << max_diff_index);
        //     INFO("Old value: " << std::setprecision(cout_precision)
        //                        << old_out_array(max_diff_index, ivar)
        //                        << ", Current value: "
        //                        << current_out_array(max_diff_index, ivar));
        //     CHECK(max_diff == doctest::Approx(0.0).epsilon(1e-12));
        // }

        // // GPU barrier
        // amrex::Gpu::streamSynchronize();
    }
    amrex::Finalize();
}
