/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// Doctest header
#include "doctest.h"

// Test include
#include "LagrangeUnitTest.hpp"

// Common includes
#include "doctestCLIArgs.hpp"

// AMReX includes
#include "AMReX.H"
#include <AMReX_Array.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Particles.H>
#include <AMReX_Print.H>
#include "AMReX_IntVect.H"

// Base includes
#include "DefaultLevelFactory.hpp"
#include "GRAMR.hpp"
#include "GRParmParse.hpp"

// Problem specific includes
#include "Derivative.hpp"
#include "DerivativeSetup.hpp"
#include "InterpolatorTestLevel.hpp"
#include "ParticleInterpolators.hpp"
#include "PolynomialTest.hpp"
#include "SimulationParameters.hpp"

// Others
#include <filesystem>

// An interpolation problem borrowed from original GRChombo tests (using only
// one polynomial example however) We treat the polynomial as a derived
// variable, and then interpolate it to some points using the
// ParticleInterpolators class.

void run_lagrange_test()
{
    // Use an input file that is in the same directory as this file for the
    // second argument
    std::filesystem::path this_file(__FILE__);
    std::filesystem::path input_file =
        this_file.parent_path() /
        std::filesystem::path("AMRInterpolatorTest.inputs");
    char *input_file_c_str = strdup(input_file.c_str());

    auto new_args = doctest::cli_args;
    new_args.insert(1, input_file_c_str);

    int new_argc    = new_args.argc();
    char **new_argv = new_args.argv();

    // NOLINTNEXTLINE(bugprone-casting-through-void) // Open MPI triggers this
    amrex::Initialize(new_argc, new_argv);
    {
        // Simulation parameters
        GRParmParse pp;
        SimulationParameters sim_params(pp);
        GRAMR::set_simulation_parameters(sim_params);

        // Set the center
        PolynomialTest::set_center(sim_params.center);

        // Set up the AMR object
        DefaultLevelFactory<InterpolatorTestLevel> interpolator_test_level_fact;
        GRAMR gr_amr(&interpolator_test_level_fact);
        gr_amr.init(0., sim_params.stop_time);

        // Build the polynomial on the grid; we iterate over levels and fill a
        // MultiFab at each level
        const int finest = gr_amr.finestLevel(); // get finest level here
        const int ngrow  = 2;                    // no of ghost cells
        std::vector<std::unique_ptr<amrex::MultiFab>> poly_by_lev(finest + 1);

        // iterate
        for (int lev = 0; lev <= finest; ++lev)
        {
            auto &L        = gr_amr.getLevel(lev);       // level
            auto &state    = L.get_new_data(State_Type); // state data
            const auto &ba = state.boxArray();           // box array
            const auto &dm = state.DistributionMap();    // distribution map

            poly_by_lev[lev] = std::make_unique<amrex::MultiFab>(
                ba, dm, 1, ngrow); // we have 1 number of comps

            // fill the MultiFab: note that destination comp = 0, number of
            // components = 1, geometry is read from AMRLevel see e.g. here
            // (https://amrex-codes.github.io/amrex/doxygen/classamrex_1_1AmrLevel.html)
            PolynomialTest::compute_mf(
                *poly_by_lev[lev], 0, 1, *poly_by_lev[lev], // unused
                L.Geom(), /*time=*/0.0, /*bcrec=*/nullptr, lev);
        }
        amrex::Gpu::streamSynchronize();

        // Build the point from sim_params
        const int num_points = sim_params.num_points;

        std::vector<double> A(num_points);
        std::vector<double> interp_x(num_points);
        std::vector<double> interp_y(num_points);
        std::vector<double> interp_z(num_points);

        double extract_radius = sim_params.L / 4;

        for (int ipoint = 0; ipoint < num_points; ++ipoint)
        {
            double phi   = ipoint * 2. * M_PI / num_points;
            double theta = ipoint * M_PI / num_points;
            interp_x[ipoint] =
                sim_params.center[0] + extract_radius * cos(phi) * sin(theta);
            interp_y[ipoint] =
                sim_params.center[1] + extract_radius * sin(phi) * sin(theta);
            interp_z[ipoint] =
                sim_params.center[2] + extract_radius * cos(theta);
        }

        // std::cout << "Interp_x[0] " << interp_x[0] - sim_params.center[0] <<
        // std::endl; std::cout << "Interp_y[0] " << interp_y[0] -
        // sim_params.center[1] << std::endl; std::cout << "Interp_z[0] " <<
        // interp_z[0] - sim_params.center[2] << std::endl;

        // set-up query
        InterpolationQueryParticle query(num_points);
        query.setCoords(0, interp_x.data())
            .setCoords(1, interp_y.data())
            .setCoords(2, interp_z.data())
            .addComp(0, A.data(), Derivative::LOCAL, VariableType::derived);

        // set up interpolation using Particles
        ParticleInterpolators<1> interpolator(sim_params.boundary_params, 0);
        interpolator.set_gramr_ptr(&gr_amr);
        interpolator.populate_from_query(query);

        std::vector<const amrex::MultiFab *> fields;
        fields.reserve(poly_by_lev.size());
        for (const auto &mf_uptr : poly_by_lev)
        {
            fields.push_back(
                mf_uptr.get()); // convert to const amrex::MultiFab* to feed to
                                // interp call like below
        }

        std::vector<int> comps(fields.size(),
                               0); // we have only one comp, which is 0

        interpolator.interpolate_to_particle_from_derived_fields(fields);

        interpolator.interp(query, VariableType::derived);

        if (amrex::ParallelDescriptor::MyProc() == 0)
        {
            double diff = 0; // absolute error

            for (int ipoint = 0; ipoint < num_points; ++ipoint)
            {
                double x = interp_x[ipoint] - sim_params.center[0];
                double y = interp_y[ipoint] - sim_params.center[1];
                double z = interp_z[ipoint] - sim_params.center[2];

                double value_A = 42. + x * x + y * y * z * z;

                diff = fabs(A[ipoint] - value_A);

                amrex::Print() << "Absolute error is " << std::setprecision(10)
                               << diff << "\n";

                CHECK(diff == doctest::Approx(0.0).epsilon(1e-10));
            }
        }
    }

    amrex::Finalize();
}
