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
#include "Cell.hpp"
#include "HarmonicTest.hpp"

enum
{
    c_phi,
    NUM_SPHERICAL_HARMONICS_VARS
};

TEST_CASE("Spherical Harmonic")
{
    amrex::Initialize(MPI_COMM_WORLD);

    const int N_GRID = 64;
    amrex::Box box(amrex::IntVect::TheZeroVector(),
                   amrex::IntVect(N_GRID - 1, N_GRID - 1, N_GRID - 1));
    amrex::FArrayBox in_fab(box, NUM_SPHERICAL_HARMONICS_VARS,
                            amrex::The_Managed_Arena());
    amrex::FArrayBox out_fab(box, NUM_SPHERICAL_HARMONICS_VARS,
                             amrex::The_Managed_Arena());
    amrex::FArrayBox diff_fab(box, NUM_SPHERICAL_HARMONICS_VARS,
                              amrex::The_Managed_Arena());
    double length = 64.0;

    const double dx     = length / (N_GRID);
    const double center = 0.5 * length;
    auto in_array       = in_fab.array();
    auto out_array      = out_fab.array();
    auto diff_array     = diff_fab.array();

    std::array<double, AMREX_SPACEDIM> center_vector = {center, center, center};
    HarmonicTest harmonic_test(center_vector, dx);

    amrex::ParallelFor(
        box,
        [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            const double x  = (i + 0.5) * dx - center;
            const double y  = (j + 0.5) * dx - center;
            const double z  = (k + 0.5) * dx - center;
            const double r  = std::max(1e-6, std::sqrt(x * x + y * y + z * z));
            const double rr = r * r;
            const double rr_inv = 1.0 / rr;
            const double rho    = std::max(1e-6, std::sqrt(x * x + y * y));

            const amrex::IntVect iv{i, j, k};
            // here testing the es = -1, el = 2, em = -1 case
            // and also the calculation of r in coords
            double harmonic = sqrt(5.0 / 16.0 / M_PI) * x *
                              (2 * z * z - z * r - rr) * rr_inv / rho;
            in_array(iv, c_phi) = harmonic * rr_inv;

            amrex::CellData<amrex::Real> cell = out_array.cellData(i, j, k);
            harmonic_test.compute(i, j, k, cell);

            diff_array(iv, c_phi) =
                std::fabs(in_array(iv, c_phi) - out_array(iv, c_phi));
        });

    amrex::Gpu::streamSynchronize();

    const int cout_precision    = Catch::StringMaker<amrex::Real>::precision;
    const double test_tolerance = 1e-14;

    amrex::Real max_diff = 0.0;
    amrex::IntVect max_diff_index{};

    diff_fab.maxIndex(box, max_diff, max_diff_index, c_phi);

    INFO("Max diff = " << std::setprecision(cout_precision) << max_diff
                       << " at " << max_diff_index);
    INFO("SphericalHarmonics computed value = "
         << std::setprecision(cout_precision)
         << out_array(max_diff_index, c_phi));
    INFO("Correct value = " << std::setprecision(cout_precision)
                            << in_array(max_diff_index, c_phi));
    CHECK_THAT(max_diff, Catch::Matchers::WithinAbs(0.0, test_tolerance));

    amrex::Finalize();
}
