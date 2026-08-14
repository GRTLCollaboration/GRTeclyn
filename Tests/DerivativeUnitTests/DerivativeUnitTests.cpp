/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// Doctest header
#include "doctest.h"

// Test header
#include "DerivativeUnitTests.hpp"

// Common includes
#include "doctestCLIArgs.hpp"

// AMReX includes
#include "AMReX.H"
#include "AMReX_FArrayBox.H"

// Other includes
#include <iostream>

// Our includes
#include "DerivativeTestsCompute.hpp"
#include "FourthOrderDerivatives.hpp"
#include "SixthOrderDerivatives.hpp"

// Helper function to validate derivative results for different orders
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
void validate_derivatives(const amrex::Box &box,
                          const amrex::FArrayBox &out_fab, const amrex::Real dx,
                          const char *order_name, amrex::Real diss_factor,
                          int diss_power)
{
    constexpr amrex::Real test_threshold = 1e-10;
    const auto &out_c_array              = out_fab.const_array();

    amrex::LoopOnCpu(
        box,
        // NOLINTNEXTLINE(readability-function-cognitive-complexity)
        [=](int ix, int iy, int iz)
        {
            // only 1 cell in the y direction
            const amrex::Real x = (0.5 + ix) * dx;
            const amrex::Real z = (0.5 + iz) * dx;

            amrex::IntVect iv(ix, iy, iz);
            const auto &cell_data = out_c_array.cellData(ix, iy, iz);

            // First order derivative tests
            CHECK_MESSAGE(
                cell_data[c_d1] ==
                    doctest::Approx(2. * x * (z - 0.5)).epsilon(test_threshold),
                "Failed diff1 scalar (", order_name, ") at ", iv);

            CHECK_MESSAGE(
                cell_data[c_d1_v3] ==
                    doctest::Approx(2. * x * (z - 0.5)).epsilon(test_threshold),
                "Failed diff1 vector (", order_name, ") at ", iv);

            CHECK_MESSAGE(
                cell_data[c_d1_t33] ==
                    doctest::Approx(2. * x * (z - 0.5)).epsilon(test_threshold),
                "Failed diff1 tensor (", order_name, ") at ", iv);

            // Second order derivative tests
            CHECK_MESSAGE(cell_data[c_d2] ==
                              doctest::Approx(2. * x).epsilon(test_threshold),
                          "Failed diff2 scalar (", order_name, ") at ", iv);

            CHECK_MESSAGE(cell_data[c_d2_v3] ==
                              doctest::Approx(2. * x).epsilon(test_threshold),
                          "Failed diff2 vector (", order_name, ") at ", iv);

            CHECK_MESSAGE(
                cell_data[c_d2_t31] ==
                    doctest::Approx(2. * (z - 0.5)).epsilon(test_threshold),
                "Failed diff2 tensor (", order_name, ") at ", iv);

            CHECK_MESSAGE(cell_data[c_d2_sym_t33] ==
                              doctest::Approx(2. * x).epsilon(test_threshold),
                          "Failed diff2 symmetric tensor (", order_name,
                          ") at ", iv);

            CHECK_MESSAGE(
                cell_data[c_d2_mixed] ==
                    doctest::Approx(2. * (z - 0.5)).epsilon(test_threshold),
                "Failed mixed diff2 (", order_name, ") at ", iv);

            // Advection and dissipation tests
            // diss_factor and diss_power are different between Fourth and Sixth
            // Order Derivatives The SixthOrderDerivatives use 8th order
            // dissipation
            amrex::Real expected_diss =
                diss_factor * (1. + z * (z - 1.)) * pow(dx, diss_power);

            CHECK_MESSAGE(
                cell_data[c_diss] ==
                    doctest::Approx(expected_diss).epsilon(test_threshold),
                "Failed dissipation (", order_name, ") at ", iv);

            CHECK_MESSAGE(
                cell_data[c_advec_down] ==
                    doctest::Approx(-2. * z * (z - 1.) - 3. * x * (2. * z - 1.))
                        .epsilon(test_threshold),
                "Failed advection down (", order_name, ") at ", iv);

            CHECK_MESSAGE(
                cell_data[c_advec_up] ==
                    doctest::Approx(2. * z * (z - 1.) + 3. * x * (2. * z - 1.))
                        .epsilon(test_threshold),
                "Failed advection up (", order_name, ") at ", iv);
        });
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
void run_derivative_unit_tests()
{
    int amrex_argc    = doctest::cli_args.argc();
    char **amrex_argv = doctest::cli_args.argv();
    // NOLINTNEXTLINE(bugprone-casting-through-void) // Open MPI triggers this
    amrex::Initialize(amrex_argc, amrex_argv);
    {
        constexpr int num_cells  = 32;
        constexpr int num_ghosts = 4;
        // box is flat in y direction to make test cheaper
        amrex::IntVect domain_hi_vect(num_cells - 1, 0, num_cells - 1);
        amrex::Box box(amrex::IntVect::TheZeroVector(), domain_hi_vect);
        amrex::Box ghosted_box = box;
        ghosted_box.grow(num_ghosts);

        amrex::FArrayBox in_fab(ghosted_box, NUM_DERIVATIVES_VARS,
                                amrex::The_Managed_Arena());
        amrex::FArrayBox out_fab(box, NUM_DERIVATIVES_VARS,
                                 amrex::The_Managed_Arena());

        const amrex::Real dx = 1.0 / num_cells;

        const amrex::Array4<amrex::Real> &in_array = in_fab.array();

        amrex::ParallelFor(ghosted_box,
                           [=] AMREX_GPU_DEVICE(int ix, int iy, int iz)
                           {
                               // no point having data varying wrt y as we only
                               // 1 true cell in that dimension
                               const amrex::Real x = (0.5 + ix) * dx;
                               const amrex::Real z = (0.5 + iz) * dx;
                               for (int ivar = 0; ivar < in_array.nComp();
                                    ++ivar)
                               {
                                   in_array(ix, iy, iz, ivar) = x * z * (z - 1);
                               }
                               // The dissipation component is special:
                               in_array(ix, iy, iz, c_diss) =
                                   (pow(z - 0.5, 6) - 0.015625) / 720. +
                                   (z - 1) * z * pow(x, 6) / 720.;
                           });

        amrex::Gpu::streamSynchronize();
        AMREX_GPU_ERROR_CHECK();

        const auto &out_array  = out_fab.array();
        const auto &in_c_array = in_fab.const_array();

        SUBCASE("Fourth order derivatives")
        {
            DerivativeTestsCompute<FourthOrderDerivatives>
                derivative_tests_compute(dx);
            amrex::ParallelFor(box,
                               [=] AMREX_GPU_DEVICE(int ix, int iy, int iz)
                               {
                                   derivative_tests_compute(
                                       ix, iy, iz, out_array, in_c_array);
                               });

            amrex::Gpu::streamSynchronize();
            AMREX_GPU_ERROR_CHECK();

            validate_derivatives(box, out_fab, dx, "fourth order", 1.0 / 64.0,
                                 5);
        }

        SUBCASE("Sixth order derivatives")
        {
            DerivativeTestsCompute<SixthOrderDerivatives>
                derivative_tests_compute(dx);
            amrex::ParallelFor(box,
                               [=] AMREX_GPU_DEVICE(int ix, int iy, int iz)
                               {
                                   derivative_tests_compute(
                                       ix, iy, iz, out_array, in_c_array);
                               });

            amrex::Gpu::streamSynchronize();
            AMREX_GPU_ERROR_CHECK();

            validate_derivatives(box, out_fab, dx, "sixth order", -1.0 / 256.0,
                                 7);
        }
    }
    amrex::Finalize();
}