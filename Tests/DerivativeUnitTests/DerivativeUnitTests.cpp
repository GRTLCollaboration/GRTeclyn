/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// Catch2 header
#include "catch_amalgamated.hpp"

// AMReX includes
#include "AMReX.H"
#include "AMReX_FArrayBox.H"

// Other includes
#include <iostream>

// Our includes
#include "DerivativeTestsCompute.hpp"
#include "FourthOrderDerivatives.hpp"
// #include "SixthOrderDerivatives.hpp"

TEST_CASE("Derivatives")
{
    amrex::Initialize(MPI_COMM_WORLD);
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

        const double dx = 1.0 / num_cells;

        const amrex::Array4<amrex::Real> &in_array = in_fab.array();

        amrex::ParallelFor(ghosted_box,
                           [=] AMREX_GPU_DEVICE(int i, int j, int k)
                           {
                               // no point having data varying wrt y as we only
                               // 1 true cell in that dimension
                               const double x = (0.5 + i) * dx;
                               const double z = (0.5 + k) * dx;
                               for (int ivar = 0; ivar < in_array.nComp();
                                    ++ivar)
                               {
                                   in_array(i, j, k, ivar) = x * z * (z - 1);
                               }
                               // The dissipation component is special:
                               in_array(i, j, k, c_diss) =
                                   (pow(z - 0.5, 6) - 0.015625) / 720. +
                                   (z - 1) * z * pow(x, 6) / 720.;
                           });

        amrex::Gpu::streamSynchronize();

        const auto &out_array  = out_fab.array();
        const auto &in_c_array = in_fab.const_array();

        SECTION("Fourth order derivatives")
        {
            DerivativeTestsCompute<FourthOrderDerivatives>
                derivative_tests_compute(dx);
            amrex::ParallelFor(
                box, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                { derivative_tests_compute(i, j, k, out_array, in_c_array); });

            amrex::Gpu::streamSynchronize();

            constexpr amrex::Real test_threshold = 1e-10;

            const auto &out_c_array = out_fab.const_array();

            amrex::LoopOnCpu(
                box,
                [=](int i, int j, int k)
                {
                    // only 1 cell in the y direction
                    const double x = (0.5 + i) * dx;
                    const double z = (0.5 + k) * dx;

                    amrex::IntVect iv(i, j, k);

                    DerivativeTestsCompute<FourthOrderDerivatives>::Vars<
                        amrex::Real>
                        vars;
                    const auto &cell_data = out_c_array.cellData(i, j, k);
                    load_vars(cell_data, vars);

                    INFO("diff1 (fourth order) at " << iv);
                    CHECK_THAT(vars.d1,
                               Catch::Matchers::WithinAbs(2. * x * (z - 0.5),
                                                          test_threshold));

                    INFO("diff2 (fourth order) at " << iv);
                    CHECK_THAT(vars.d2, Catch::Matchers::WithinAbs(
                                            2. * x, test_threshold));

                    INFO("mixed diff2 (fourth order) at " << iv);
                    CHECK_THAT(vars.d2_mixed,
                               Catch::Matchers::WithinAbs(2. * (z - 0.5),
                                                          test_threshold));

                    INFO("dissipation (fourth order) at " << iv);
                    CHECK_THAT(vars.diss,
                               Catch::Matchers::WithinAbs((1. + z * (z - 1.)) *
                                                              pow(dx, 5) / 64.,
                                                          test_threshold));

                    INFO("advection down (fourth order) at " << iv);
                    CHECK_THAT(vars.advec_down,
                               Catch::Matchers::WithinAbs(
                                   -2. * z * (z - 1.) - 3. * x * (2. * z - 1.),
                                   test_threshold));

                    INFO("advection up (fourth order) at " << iv);
                    CHECK_THAT(vars.advec_up,
                               Catch::Matchers::WithinAbs(
                                   2. * z * (z - 1.) + 3. * x * (2. * z - 1.),
                                   test_threshold));
                });
        }

        // SECTION("Sixth order derivatives")
        // {
        //     DerivativeTestsCompute<SixthOrderDerivatives>
        //     derivative_tests_compute(
        //         dx);
        //     amrex::ParallelFor(
        //         box, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        //         { derivative_tests_compute(i, j, k, out_array, in_c_array);
        //         });

        //     amrex::Gpu::streamSynchronize();

        //     constexpr amrex::Real test_threshold = 1e-10;

        //     const auto &out_c_array = out_fab.const_array();

        //     amrex::LoopOnCpu(
        //         box,
        //         [=] (int i, int j, int k)
        //         {
        //             // only 1 cell in the y direction
        //             const double x = (0.5 + i) * dx;
        //             const double z = (0.5 + k) * dx;

        //             amrex::IntVect iv(i, j, k);

        //             DerivativeTestsCompute<SixthOrderDerivatives>::Vars<
        //                 amrex::Real>
        //                 vars;
        //             const auto &cell_data = out_c_array.cellData(i, j, k);
        //             load_vars(cell_data, vars);

        //             INFO("diff1 (sixth order) at " << iv);
        //             CHECK_THAT(vars.d1, Catch::Matchers::WithinAbs(
        //                                     2. * x * (z - 0.5),
        //                                     test_threshold));

        //             INFO("diff2 (sixth order) at " << iv);
        //             CHECK_THAT(vars.d2,
        //                        Catch::Matchers::WithinAbs(2. * x,
        //                        test_threshold));

        //             INFO("mixed diff2 (sixth order) at " << iv);
        //             CHECK_THAT(vars.d2_mixed, Catch::Matchers::WithinAbs(
        //                                           2. * (z - 0.5),
        //                                           test_threshold));

        //             INFO("dissipation (sixth order) at " << iv);
        //             CHECK_THAT(vars.diss,
        //                        Catch::Matchers::WithinAbs((1. + z * (z - 1.))
        //                        *
        //                                                       pow(dx, 5)
        //                                                       / 64.,
        //                                                   test_threshold));

        //             INFO("advection down (sixth order) at " << iv);
        //             CHECK_THAT(vars.advec_down,
        //                        Catch::Matchers::WithinAbs(
        //                            -2. * z * (z - 1.) - 3. * x * (2. * z
        //                            - 1.), test_threshold));

        //             INFO("advection up (sixth order) at " << iv);
        //             CHECK_THAT(vars.advec_up,
        //                        Catch::Matchers::WithinAbs(
        //                            2. * z * (z - 1.) + 3. * x * (2. * z
        //                            - 1.), test_threshold));
        //         });
        // }
    }
    amrex::Finalize();
}