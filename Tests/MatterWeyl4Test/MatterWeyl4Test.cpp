/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */
#define COMPARE_WITH_CHF
#define COVARIANTZ4

// AMReX includes
#include <AMReX.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_FabArray.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>

#ifdef AMREX_USE_HDF
#include "H5Cpp.h"
#include <AMReX_PlotFileUtilHDF5.H>
using namespace H5;
#endif

// Doctest includes
#include "doctest.h"
#include "doctestCLIArgs.hpp"

#include <iomanip>
#include <iostream>
#include <sys/time.h>

// common includes
#include <InitialData.hpp> //includes StateVariables.hpp

// test header
#include "MatterWeyl4Test.hpp"

#ifdef AMREX_USE_HDF5
#include "H5Cpp.h"
#include <AMReX_PlotFileUtilHDF5.H>
using namespace H5;
#endif

// GRTeclyn includes

#include "CCZ4RHS.hpp"
#include "DefaultPotential.hpp"
#include "EMTensor.hpp"
#include "MatterCCZ4RHS.hpp"
#include "MatterWeyl4.hpp"
#include "ScalarField.hpp"
// #include "Weyl4.hpp"
#include "simd.hpp"
#include <array>

void run_matter_weyl4_test()
{
    int amrex_argc    = doctest::cli_args.argc();
    char **amrex_argv = doctest::cli_args.argv();

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
        amrex::DistributionMapping distribution_mapping{box_array};
        amrex::MFInfo mf_info;
        mf_info.SetArena(amrex::The_Managed_Arena());

        amrex::MultiFab in_mf{box_array, distribution_mapping, NUM_VARS,
                              num_ghosts, mf_info};

        amrex::FArrayBox in_fab{ghosted_box, NUM_VARS,
                                amrex::The_Managed_Arena()};

        const auto &in_arrays = in_mf.arrays();
        const auto &in_array  = in_fab.array();

        amrex::ParallelFor(
            ghosted_box,
            [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                const amrex::IntVect iv{i, j, k};
                const amrex::RealVect coords = amrex::RealVect{iv} * dx;
                amrex::Real x                = coords[0];
                amrex::Real y                = coords[1];
                amrex::Real z                = coords[2];

                random_ccz4_initial_data(iv, in_array, coords);

                // Theta is zero for BSSN
                in_array(i, j, k, c_Theta) = 0.0;

                // the initial data doesn't include phi or Pi so do it here:
                in_array(i, j, k, c_phi) = 0.21232 * sin(x * 2.1232 * 3.14) *
                                           cos(y * 2.5123 * 3.15) *
                                           cos(z * 2.1232 * 3.14);
                in_array(i, j, k, c_Pi) = 0.4112 * sin(x * 4.123 * 3.14) *
                                          cos(y * 2.2312 * 3.15) *
                                          cos(z * 2.5123 * 3.14);
            });

        amrex::Gpu::streamSynchronize();

        // Setup scalar field calculations

        using DefaultScalarField = ScalarField<DefaultPotential>;

        ScalarField<DefaultPotential> my_scalar_field(DefaultPotential());

        // set up weyl4 calculation

        constexpr int dcomp_weyl4 = 0;
        constexpr int num_comps_weyl4 =
            2; // compute will automatically +1 for imaginary component
        double G_Newton = 1.0;
        std::array<double, AMREX_SPACEDIM> center{0.0, 0.0, 0.0};

        MatterWeyl4<DefaultScalarField> matter_weyl4(
            DefaultScalarField(DefaultPotential()), center, dx, dcomp_weyl4,
            CCZ4RHS<>::USE_BSSN, G_Newton);

        // Constructor for EMTensor
        constexpr int dcomp_rho = num_comps_weyl4;
        EMTensor<DefaultScalarField> scalar_field_emtensor(
            DefaultScalarField(DefaultPotential()), dx, dcomp_rho);

        constexpr int num_comps =
            dcomp_rho + 1; // just Weyl4_Re, Weyl4_Im, rho, don't bother
        //            storing Si, Sij

        amrex::MultiFab out_mf{box_array, distribution_mapping, num_comps,
                               num_ghosts, mf_info};

        amrex::FArrayBox out_fab{box, num_comps, amrex::The_Managed_Arena()};

        const auto &in_c_arrays = in_mf.const_arrays();
        const auto &out_arrays  = out_mf.arrays();

        const auto &in_c_array = in_fab.const_array();
        const auto &out_array  = out_fab.array();

        amrex::ParallelFor(
            box,
            [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                matter_weyl4.compute(i, j, k, out_array, in_c_array);

                // Rho is also computed here
                scalar_field_emtensor.compute(i, j, k, out_array, in_c_array);
            });

#if AMREX_USE_HDF5
        int coord_sys                        = 0;
        amrex::Vector<std::string> var_names = {"Weyl4_Re", "Weyl4_Im", "rho"};

        const H5std_string grteclyn_hdf5_file =
            "MatterWeyl4Test/MatterWeyl4TestOut.h5";

        // open the hdf5 file for writing
        hid_t fid =
            H5Fcreate(grteclyn_hdf5_file.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                      H5P_DEFAULT); // H5F_ACC_TRUNC = if file exists open with
                                    // read/write access, otherwise create file

        // // create the group
        char level_name[8] = "level_0"; // only the one level

        hid_t grp =
            H5Gcreate(fid, level_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        for (int i = 0; i < num_comps; i++)
        {

            // // create the dataset
            char dataname[32] = "data:datatype=0";
            hsize_t hs_procsize[1]; // only 1 dimension
            hs_procsize[0] = box.length(0) * box.length(1) *
                             box.length(2); // print 1 variable at a time

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
        H5Fclose(fid);

        std::cout.flush();

        std::string h5diff_tol = "1e-10";
        std::string grchombo_weyl4_hdf5_file =
            "MatterWeyl4Test/GRChomboMatterWeyl4Test.h5";

        // this is saved separately because the GRChombo rho is outputted to the
        // input fab and these have ghost cells while we want to compare the
        // result to the GRTeclyn output fab which doesn't
        std::string grchombo_emtensor_hdf5_file =
            "MatterWeyl4Test/GRChomboEMTensorTest.h5";

        std::string hdf5_internal_path_to_rho      = "/level_0/rho";
        std::string hdf5_internal_path_to_Weyl4_Re = "/level_0/Weyl4_Re";
        std::string hdf5_internal_path_to_Weyl4_Im = "/level_0/Weyl4_Im";

        std::string h5diff_command  = "h5diff";
        h5diff_command             += " -d " + h5diff_tol;

        std::string h5diff_command_weyl4 = h5diff_command + " " +
                                           grteclyn_hdf5_file + " " +
                                           grchombo_weyl4_hdf5_file;
        std::string h5diff_command_weyl4_real =
            h5diff_command_weyl4 + " " + hdf5_internal_path_to_Weyl4_Re + " " +
            hdf5_internal_path_to_Weyl4_Re;

        std::string h5diff_command_weyl4_im =
            h5diff_command_weyl4 + " " + hdf5_internal_path_to_Weyl4_Im + " " +
            hdf5_internal_path_to_Weyl4_Im;

        std::string h5diff_command_emtensor = h5diff_command + " " +
                                              grteclyn_hdf5_file + " " +
                                              grchombo_emtensor_hdf5_file;
        h5diff_command_emtensor +=
            " " + hdf5_internal_path_to_rho + " " + hdf5_internal_path_to_rho;

        amrex::Print() << "Test command (Weyl4 - real): "
                       << h5diff_command_weyl4_real << std::endl;

        amrex::Print() << "Test command (Weyl4 - imaginary): "
                       << h5diff_command_weyl4_im << std::endl;

        amrex::Print() << "Test command (EMTensor): " << h5diff_command_emtensor
                       << std::endl;

        int h5diff_status_weyl4_re =
            std::system(h5diff_command_weyl4_real.c_str());
        int h5diff_status_weyl4_im =
            std::system(h5diff_command_weyl4_im.c_str());
        int h5diff_status_emtensor =
            std::system(h5diff_command_emtensor.c_str());

        int h5diff_retval = -1;

        // Use POSIX macros to get the exit code
        if (WIFEXITED(h5diff_status_weyl4_re))
        {
            h5diff_retval = WEXITSTATUS(h5diff_status_weyl4_re);
        }

        CHECK(h5diff_retval == 0);

        h5diff_retval = -1;

        if (WIFEXITED(h5diff_status_weyl4_im))
        {
            h5diff_retval = WEXITSTATUS(h5diff_status_weyl4_im);
        }

        CHECK(h5diff_retval == 0);

        h5diff_retval = -1;

        if (WIFEXITED(h5diff_status_emtensor))
        {
            h5diff_retval = WEXITSTATUS(h5diff_status_emtensor);
        }

        CHECK(h5diff_retval == 0);

#endif
    }
    amrex::Finalize();
}
