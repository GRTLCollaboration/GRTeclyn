/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// Doctest header
#include "doctest.h"

// AMReX includes
#include "AMReX.H"
#include "AMReX_FArrayBox.H"

// Other includes
#include "PositiveChiAndAlpha.hpp"
#include "Tensor.hpp"

TEST_CASE("Positive Chi and Alpha")
{
    amrex::Initialize(MPI_COMM_WORLD);
    {
        constexpr int N_GRID = 8;
        amrex::Box box(amrex::IntVect(0, 0, 0),
                       amrex::IntVect(N_GRID - 1, N_GRID - 1, N_GRID - 1));
        amrex::FArrayBox in_fab(box, NUM_VARS, amrex::The_Managed_Arena());

        const amrex::Array4<amrex::Real> &in_array = in_fab.array();

        amrex::ParallelFor(box,
                           [=] AMREX_GPU_DEVICE(int ix, int iy, int iz)
                           {
                               const amrex::IntVect iv{ix, iy, iz};
                               double value = (ix < N_GRID / 2) ? 1 : 1e-10;

                               in_array(iv, c_chi)   = value;
                               in_array(iv, c_lapse) = value;
                           });

        amrex::Gpu::streamSynchronize();

        amrex::ParallelFor(box,
                           [=] AMREX_GPU_DEVICE(int ix, int iy, int iz)
                           {
                               auto cell = in_array.cellData(ix, iy, iz);
                               PositiveChiAndAlpha()(cell);
                           });

        amrex::Gpu::streamSynchronize();

        constexpr double test_threshold = 1e-15;

        // We have to do this on the host as are using Catch2 functions
        amrex::LoopOnCpu(
            box,
            [=](int ix, int iy, int iz)
            {
                const amrex::IntVect iv(ix, iy, iz);

                double correct_value = (ix < N_GRID / 2) ? 1 : 1e-4;
                INFO("At " << iv);
                CHECK(in_array(iv, c_chi) ==
                      doctest::Approx(correct_value).epsilon(test_threshold));
                CHECK(in_array(iv, c_lapse) ==
                      doctest::Approx(correct_value).epsilon(test_threshold));
            });
    }
    amrex::Finalize();
}
