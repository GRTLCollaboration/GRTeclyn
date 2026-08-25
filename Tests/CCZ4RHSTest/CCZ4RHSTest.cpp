/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// Doctest header
#include "doctest.h"

// Test header
#include "CCZ4RHSTest.hpp"

// Common test headers
#include "InitialData.hpp"
#include "doctestCLIArgs.hpp"

// GRTeclyn headers
#include "CCZ4RHS.hpp"
#include "FourthOrderDerivatives.hpp"

// Old GRTeclyn headers for comparison
#include "CCZ4RHS-fdf5a7a.hpp"
#include "FourthOrderDerivatives-fdf5a7a.hpp"

// AMReX headers
#include "AMReX.H"
#include "AMReX_FArrayBox.H"

void run_ccz4_rhs_test()
{
    int amrex_argc    = doctest::cli_args.argc();
    char **amrex_argv = doctest::cli_args.argv();
    // NOLINTNEXTLINE(bugprone-casting-through-void) // Open MPI triggers this
    amrex::Initialize(amrex_argc, amrex_argv);
    {
        constexpr int num_cells  = 32;
        constexpr int num_ghosts = 3;
        constexpr amrex::Real dx = 0.5 / ((amrex::Real)num_cells - 1.0);

        amrex::Box box(
            amrex::IntVect(0, 0, 0),
            amrex::IntVect(num_cells - 1, num_cells - 1, num_cells - 1));

        amrex::Box ghosted_box = box;
        ghosted_box.grow(num_ghosts);

        amrex::FArrayBox in_fab{ghosted_box, NUM_CCZ4_VARS,
                                amrex::The_Managed_Arena()};

        const amrex::Array4<amrex::Real> &in_array = in_fab.array();
        amrex::ParallelFor(ghosted_box,
                           [=] AMREX_GPU_DEVICE(int ix, int iy, int iz)
                           {
                               const amrex::IntVect iv{ix, iy, iz};
                               const amrex::RealVect coords =
                                   amrex::RealVect{iv} * dx;

                               random_ccz4_initial_data(iv, in_array, coords);
                           });

        // These need to be const and so declared separately
        const int use_covariantZ4 = 1;
        const int formulation     = 0;

        GRParmParse pp;
        pp.add("ccz4.kappa1", 0.1);
        pp.add("ccz4.kappa2", 0.0);
        pp.add("ccz4.kappa3", 1.0);
        pp.add("ccz4.covariantZ4", use_covariantZ4);

        pp.add("gauge.shift_Gamma_coeff", 0.75);
        pp.add("gauge.lapse_advec_coeff", 0.0);
        pp.add("gauge.lapse_power", 1.0);
        pp.add("gauge.lapse_coeff", 2.0);
        pp.add("gauge.shift_advec_coeff", 0.0);
        pp.add("gauge.eta", 1.82);

        pp.add("evolution.sigma", 0.3);
        pp.add("ccz4.formulation", formulation);

        Old::CCZ4_params_t<Old::MovingPunctureGauge::params_t> old_ccz4_params;
        old_ccz4_params.kappa1            = 0.1;
        old_ccz4_params.kappa2            = 0;
        old_ccz4_params.kappa3            = 1;
        old_ccz4_params.covariantZ4       = use_covariantZ4;
        old_ccz4_params.lapse_advec_coeff = 0.0;
        old_ccz4_params.lapse_power       = 1.0;
        old_ccz4_params.lapse_coeff       = 2.0;
        old_ccz4_params.shift_Gamma_coeff = 0.75;
        old_ccz4_params.shift_advec_coeff = 0;
        old_ccz4_params.eta               = 1.82;

        amrex::Real sigma = 0.3;

        CCZ4RHS<MovingPunctureGauge, FourthOrderDerivatives> current_ccz4_rhs{
            dx};

        Old::CCZ4RHS<Old::MovingPunctureGauge, Old::FourthOrderDerivatives>
            old_ccz4_rhs{old_ccz4_params, dx, sigma};

        amrex::FArrayBox current_out_fab{box, NUM_CCZ4_VARS,
                                         amrex::The_Managed_Arena()};
        amrex::FArrayBox old_out_fab{box, NUM_CCZ4_VARS,
                                     amrex::The_Managed_Arena()};
        amrex::FArrayBox diff_fab{box, NUM_CCZ4_VARS,
                                  amrex::The_Managed_Arena()};

        const auto &in_c_array        = in_fab.const_array();
        const auto &current_out_array = current_out_fab.array();
        const auto &old_out_array     = old_out_fab.array();
        const auto &diff_array        = diff_fab.array();

        // Do the current and old CCZ4RHS calculation in the same loop

        amrex::ParallelFor(
            box, [=] AMREX_GPU_DEVICE(int ix, int iy, int iz)
            { old_ccz4_rhs.compute(ix, iy, iz, old_out_array, in_c_array); });

        // The RHS is split into three different kernels

        // NOLINTBEGIN(bugprone-easily-swappable-parameters)
        amrex::ParallelFor(box,
                           [=] AMREX_GPU_DEVICE(int ix, int iy, int iz)
                           {
                               current_ccz4_rhs.compute_chi_and_h_ij(
                                   ix, iy, iz, current_out_array, in_c_array);
                           });

        amrex::ParallelFor(
            box,
            [=] AMREX_GPU_DEVICE(int ix, int iy, int iz)
            {
                current_ccz4_rhs.compute_A_ij_and_Theta_and_Gamma<
                    formulation, use_covariantZ4>(ix, iy, iz, current_out_array,
                                                  in_c_array);
            });

        amrex::ParallelFor(box,
                           [=] AMREX_GPU_DEVICE(int ix, int iy, int iz)
                           {
                               current_ccz4_rhs.calculate_gauge_rhs(
                                   ix, iy, iz, current_out_array, in_c_array);
                               current_ccz4_rhs.apply_dissipation(
                                   ix, iy, iz, current_out_array, in_c_array);
                           });

        // GPU barrier
        amrex::Gpu::streamSynchronize();

        amrex::ParallelFor(box,
                           [=] AMREX_GPU_DEVICE(int ix, int iy, int iz)
                           {
                               for (int ivar = 0; ivar < NUM_CCZ4_VARS; ++ivar)
                               {
                                   diff_array(ix, iy, iz, ivar) = std::fabs(
                                       current_out_array(ix, iy, iz, ivar) -
                                       old_out_array(ix, iy, iz, ivar));
                               }
                           });

        // NOLINTEND(bugprone-easily-swappable-parameters)

        // GPU barrier
        amrex::Gpu::streamSynchronize();

        amrex::Real max_diff = 0.0;
        amrex::IntVect max_diff_index{};

        double test_threshold = 1e-11;

        const int cout_precision = 17;
        for (int ivar = 0; ivar < NUM_CCZ4_VARS; ++ivar)
        {
            // NOLINTNEXTLINE(bugprone-chained-comparison)
            diff_fab.maxIndex<amrex::RunOn::Device>(box, max_diff,
                                                    max_diff_index, ivar);

            INFO("Max diff for var " << StateVariables::names[ivar] << ": "
                                     << std::setprecision(cout_precision)
                                     << max_diff << " at " << max_diff_index);
            INFO("Old value: " << std::setprecision(cout_precision)
                               << old_out_array(max_diff_index, ivar)
                               << ", Current value: "
                               << current_out_array(max_diff_index, ivar));
            CHECK(max_diff == doctest::Approx(0.0).epsilon(test_threshold));
        }

        // GPU barrier
        amrex::Gpu::streamSynchronize();
    }
    amrex::Finalize();
}
