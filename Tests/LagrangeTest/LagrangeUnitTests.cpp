/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// Doctest header
#include "doctest.h"

// Test include
#include "LagrangeTest.hpp"

// Common includes
#include "doctestCLIArgs.hpp"

// AMReX includes
#include "AMReX.H"
#include <AMReX_Array.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Particles.H>
#include <AMReX_Print.H>
#include "AMReX_IntVect.H"

// Our includes
#include "LagrangeInterpolation.hpp"
#include "PolynomialTest.hpp"
// #include "LinearInterpolation.hpp"

enum
{
    c_poly,
    NUM_POLYNOMIAL_VARS
};

// A made-up interpolation problem
void run_lagrange_test()
{
    int amrex_argc    = doctest::cli_args.argc();
    char **amrex_argv = doctest::cli_args.argv();
    // NOLINTNEXTLINE(bugprone-casting-through-void) // Open MPI triggers this
    amrex::Initialize(amrex_argc, amrex_argv);
    {
        // First custom thing below
        const int nx        = 32;           // number of cells
        double length       = 32.;          // total length of the grid
        const double center = 0.5 * length; // grid center
        std::array<double, AMREX_SPACEDIM> center_vector = {
            center, center, center};       // center of the grid
        const double dx     = length / nx; // grid spacing
        const double inv_dx = 1. / dx;     // inverse grid spacing

        // Build grid data
        amrex::Box box(amrex::IntVect(0, 0, 0),
                       amrex::IntVect(nx - 1, nx - 1, nx - 1));

        // define number of components in the FArrayBox (we have only 1)
        amrex::FArrayBox out_fab(
            box, NUM_POLYNOMIAL_VARS,
            amrex::The_Managed_Arena()); // managed memory arena supports both
                                         // CPU and GPU?

        auto out_array = out_fab.array(); // get the element-wise access to the
                                          // data stored in the FArrayBox

        // Ok, now we are ready to populate the polynomial on the grid
        PolynomialTest poly_test(center_vector, dx);

        amrex::ParallelFor(box,
                           [=] AMREX_GPU_DEVICE(int i, int j, int k)
                           {
                               // We are cell-centered, so...
                               const double x = (i + 0.5) * dx - center;
                               const double y = (j + 0.5) * dx - center;
                               const double z = (k + 0.5) * dx - center;

                               const amrex::IntVect iv{i, j, k};

                               double val = poly_test.compute(i, j, k);
                               out_array(i, j, k, c_poly) = val;
                           });

        amrex::Gpu::streamSynchronize();

        // Random point where we will interpolate later on below, this is also
        // the position where we put the particle on
        amrex::IntVect cell(20, 20, 20);
        double x_interp = (cell[0] + 0.1) * dx - center;
        double y_interp = (cell[1] + 0.5) * dx - center;
        double z_interp = (cell[2] + 0.7) * dx - center;

        // Get the expected value; will be used in the check later on
        double expected_val =
            poly_test.compute_polynomial(x_interp, y_interp, z_interp);

        // Define a particle container: we have only one particle, so this is
        // redundant, but will be useful for the future
        using MyParticleContainer = amrex::ParticleContainer<AMREX_SPACEDIM, 1>;
        MyParticleContainer particles;

        // Now, we are ready to set-up our particle. This requires defining: (i)
        // geometry (ii) distribution mapping and (iii) a box array.

        // First, geometry:
        amrex::RealVect prob_lo{-0.5 * length};
        amrex::RealVect prob_hi{+0.5 * length};
        amrex::RealBox real_box{
            prob_lo.dataPtr(),
            prob_hi.dataPtr()}; // real box is essentially the physical box,
                                // unlike just the 'box' in Amrex.
        int coord_sys = 0;      // for Cartesian coordinates
        amrex::Geometry geom{box, &real_box, coord_sys};

        // Second: box array
        amrex::BoxArray box_array{box};

        // Third: distribution mapping
        amrex::DistributionMapping distribution_mapping(box_array);
        particles.Define(geom, distribution_mapping, box_array);

        // For the particle:
        int lev  = 0;
        int grid = 0;
        int tile = 0;

        auto &particle_tile =
            particles.DefineAndReturnParticleTile(lev, grid, tile);
        particle_tile.resize(1); // create space for one particle
        auto &arr_particles = particle_tile.GetArrayOfStructs();
        auto &p             = arr_particles[0];

        p.id()   = MyParticleContainer::ParticleType::NextID();
        p.cpu()  = amrex::ParallelDescriptor::MyProc();
        p.pos(0) = x_interp;
        p.pos(1) = y_interp;
        p.pos(2) = z_interp;

        // For debugging
        // std::cout << "Particle initialized with ID: " << p.id() << "\n";
        // std::cout << "Particle position: (" << p.pos(0) << ", "
        //           << p.pos(1) << ", " << p.pos(2) << ")\n";

        // Now we can interpolate using 4th order Lagrange
        amrex::GpuArray<amrex::Real, 3> plo = {-0.5 * length, -0.5 * length,
                                               -0.5 * length};
        amrex::GpuArray<amrex::Real, 3> dxi = {inv_dx, inv_dx, inv_dx};
        amrex::IntVect is_nodal{0, 0, 0};

        LagrangeInterpolator<5> interp;
        // LinearInterpolator interp;
        interp.compute_weights(p, plo, dxi, is_nodal);
        amrex::ParticleReal result[1]; // One component
        amrex::Array4<const double> const_out_array = out_array;
        interp.interpolate(&const_out_array, result, c_poly, 1);

        p.rdata(0) = result[0];

        double interp_val = p.rdata(0);
        double error      = std::abs(interp_val - expected_val);

        const double test_tolerance = 1e-10; // set your desired tolerance
        const int cout_precision    = 17;

        amrex::Print() << "Interpolated value is "
                       << std::setprecision(cout_precision) << interp_val
                       << "\n"
                       << "Expected value is "
                       << std::setprecision(cout_precision) << expected_val
                       << "\n"
                       << "Absolute error is "
                       << std::setprecision(cout_precision) << error << "\n";

        CHECK(error == doctest::Approx(0.0).epsilon(test_tolerance));

        amrex::Gpu::streamSynchronize();
    }
    amrex::Finalize();
}
