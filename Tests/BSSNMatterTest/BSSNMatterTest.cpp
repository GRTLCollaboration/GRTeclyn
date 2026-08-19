/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// Doctest header
#include "doctest.h"

// Test header
#include "BSSNMatterTest.hpp"

// Common test headers
#include "InitialData.hpp"
#include "doctestCLIArgs.hpp"

// GRTeclyn headers
#include "CCZ4RHSWithMatter.hpp"
#include "ConstraintsWithMatter.hpp"
#include "DefaultPotential.hpp"
#include "GRParmParse.hpp"
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

void run_bssn_matter_test()
{
    int amrex_argc    = doctest::cli_args.argc();
    char **amrex_argv = doctest::cli_args.argv();

    // NOLINTNEXTLINE(bugprone-casting-through-void) // Open MPI triggers this
    amrex::Initialize(amrex_argc, amrex_argv);
    {

        constexpr int num_cells  = 32;
        constexpr int num_ghosts = 3;
        constexpr amrex::Real dx = 0.5 / (num_cells - 1);

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
        in_mf.setVal(0.0); // initialise to zero

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

        GRParmParse pp;
        pp.add("ccz4.kappa1", 0.0);
        pp.add("ccz4.kappa2", 0.0);
        pp.add("ccz4.kappa3", 0.0);
        pp.add("ccz4.covariantZ4", true);

        pp.add("gauge.shift_Gamma_coeff", 0.75);
        pp.add("gauge.lapse_advec_coeff", 1.0);
        pp.add("gauge.lapse_power", 1.0);
        pp.add("gauge.lapse_coeff", 2.0);
        pp.add("gauge.shift_advec_coeff", 0.0);
        pp.add("gauge.eta", 1.0);

        pp.add("evolution.sigma", 0.1);
        pp.add("ccz4.formulation", CCZ4RHS<>::USE_BSSN);

        using DefaultScalarField =
            ScalarField<DefaultPotential, FourthOrderDerivatives>;

        CCZ4RHSWithMatter<DefaultScalarField, MovingPunctureGaugeWithMatter,
                          FourthOrderDerivatives>
            current_ccz4_rhs{dx};
        MovingPunctureGaugeWithMatter moving_puncture_gauge(dx);

        // Set up the constraints
        constexpr int num_bssn_matter_vars = c_Pi + 1;
        constexpr int dcomp                = num_bssn_matter_vars;

        int num_comp_constraints = 1 + AMREX_SPACEDIM; // ham + moms
        int num_comp             = num_bssn_matter_vars + num_comp_constraints;

        amrex::MultiFab out_mf{box_array, distribution_mapping, num_comp, 0,
                               mf_info};

        amrex::FArrayBox out_fab{box, num_comp, amrex::The_Managed_Arena()};

        const auto &in_c_array    = in_mf.const_arrays();
        const auto &out_mf_array  = out_mf.arrays();
        const auto &out_fab_array = out_fab.array();

        // calculate the vacuum solution

        // NOLINTBEGIN(bugprone-easily-swappable-parameters)
        amrex::ParallelFor(
            out_mf,
            [=] AMREX_GPU_DEVICE(int ibox, int ix, int iy, int iz)
            {
                current_ccz4_rhs.compute_chi_and_h_ij(
                    ix, iy, iz, out_mf_array[ibox], in_c_array[ibox]);
            });

        amrex::ParallelFor(
            out_mf,
            [=] AMREX_GPU_DEVICE(int ibox, int ix, int iy, int iz)
            {
                current_ccz4_rhs.compute_A_ij_and_Theta_and_Gamma(
                    ix, iy, iz, out_mf_array[ibox], in_c_array[ibox]);
            });
        amrex::ParallelFor(
            out_mf,
            [=] AMREX_GPU_DEVICE(int ibox, int ix, int iy, int iz)
            {
                moving_puncture_gauge.calculate_gauge_rhs(
                    ix, iy, iz, out_mf_array[ibox], in_c_array[ibox]);
            });

        // calculate the matter contribution
        amrex::ParallelFor(
            out_mf,
            [=] AMREX_GPU_DEVICE(int ibox, int ix, int iy, int iz)
            {
                current_ccz4_rhs.add_emtensor_rhs(
                    ix, iy, iz, out_mf_array[ibox], in_c_array[ibox]);
                current_ccz4_rhs.add_matter_rhs(ix, iy, iz, out_mf_array[ibox],
                                                in_c_array[ibox]);
                current_ccz4_rhs.apply_dissipation(
                    ix, iy, iz, out_mf_array[ibox], in_c_array[ibox]);
            });

        // NOLINTEND(bugprone-easily-swappable-parameters)
        amrex::Real time = 0.0;
        int *bcrec       = nullptr;
        int level        = 0;

        ConstraintsWithMatter<DefaultScalarField>::compute_mf(
            out_mf, dcomp, num_comp_constraints, in_mf, geom, time, bcrec,
            level);

        // GPU barrier
        amrex::Gpu::streamSynchronize();

        CHECK(!out_mf.contains_nan());

#if AMREX_USE_HDF5

        amrex::Vector<std::string> bssn_matter_names(
            StateVariables::names.begin(),
            StateVariables::names.begin() + num_bssn_matter_vars);
        amrex::Vector<std::string> var_names =
            ArrayTools::concatenate(bssn_matter_names, Constraints::var_names);

        std::string grteclyn_hdf5_file = "BSSNMatterTest/BSSNMatterTest";

        // open the hdf5 file for writing
        amrex::WriteSingleLevelPlotfileHDF5(grteclyn_hdf5_file, out_mf,
                                            var_names, geom, 0.0, 0);

        std::cout.flush();

        std::string h5diff_tol         = "1e-10";
        std::string grchombo_hdf5_file = "BSSNMatterTest/"
                                         "GRChomboBSSNMatterTest.h5";
        // the GRChombo comparison file is created in
        // GRChombo/Tests/MatterCCZ4Test

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
