/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// Doctest header
#include "doctest.h"

// State vars
#include "StateVariables.hpp"

// Test include
#include "ParticleInterpolatorUnitTest.hpp"

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
#include "ParticleInterpolator.hpp"
#include "ParticleInterpolatorLevel.hpp"
#include "PolynomialDerivedQuantity.hpp"
#include "SimulationParameters.hpp"

// Others
#include <filesystem>

// An interpolation problem borrowed from original GRChombo tests. We treat one
// of the polynomials as a derived variable, and the other as a state variable.
// We interpolate them to some specified (x,y,z) points using the
// ParticleInterpolator class.

void run_particle_interpolator_test()
{
    // Use an input file that is in the same directory as this file for the
    // second argument
    std::filesystem::path this_file(__FILE__);
    std::filesystem::path input_file =
        this_file.parent_path() / std::filesystem::path("params_test.txt");
    char *input_file_c_str = strdup(input_file.c_str());

    auto new_args = doctest::cli_args;
    new_args.insert(1, input_file_c_str);

    int new_argc    = new_args.argc();
    char **new_argv = new_args.argv();

    // NOLINTNEXTLINE(bugprone-casting-through-void) // Open MPI triggers this
    amrex::Initialize(
        new_argc, new_argv,
        std::function<void()>(SimulationParameters::check_params));
    {
        // Simulation parameters
        GRParmParse pp;
        ParticleInterpolatorLevel::variableSetUp();

        // Set the center
        std::array<double, AMREX_SPACEDIM> center{};
        pp.get("geometry.center", center);
        PolynomialDerivedQuantity::set_center(center);

        // Set up the AMR object
        DefaultLevelFactory<ParticleInterpolatorLevel>
            interpolator_test_level_fact;
        GRAMR gr_amr(&interpolator_test_level_fact);

        double stop_time{};
        pp.get("evolution.stop_time", stop_time);
        gr_amr.init(0., stop_time);

        // Read from params
        int num_points{};
        pp.get("test.num_points", num_points);

        bool verbosity{};
        pp.get("particle_interpolator.verbosity", verbosity);

        std::array<double, AMREX_SPACEDIM> prob_extent{};
        pp.get("geometry.prob_extent", prob_extent);

        // Using lenght of x direction to define extraction radius
        double extract_radius = prob_extent[0] / 4;

        // Number of processes and local processes
        const int nprocs = amrex::ParallelDescriptor::NProcs();
        const int myproc = amrex::ParallelDescriptor::MyProc();

        // Partition the points across ranks
        const int base      = num_points / nprocs;
        const int remainder = num_points % nprocs;

        const int n_local =
            base + (myproc < remainder
                        ? 1
                        : 0); // local number of particles on the rank
        if (verbosity)
        {
            amrex::AllPrint() << "I am rank " << myproc << " and I query "
                              << n_local << " LOCAL particles \n";
        }
        const int start =
            myproc * base + std::min(myproc, remainder); // global start index

        // Allocate vectors for writing
        std::vector<double> A_local(n_local); // for storing derived polynomial
        std::vector<double> B_local(n_local); // for storing state polynomial
        std::vector<double> interp_x_local(n_local);
        std::vector<double> interp_y_local(n_local);
        std::vector<double> interp_z_local(n_local);

        for (int j = 0; j < n_local; ++j)
        {
            int ipoint   = start + j; // global index
            double phi   = ipoint * 2. * M_PI / num_points;
            double theta = ipoint * M_PI / num_points;

            interp_x_local[j] =
                center[0] + extract_radius * cos(phi) * sin(theta);
            interp_y_local[j] =
                center[1] + extract_radius * sin(phi) * sin(theta);
            interp_z_local[j] = center[2] + extract_radius * cos(theta);
        }

        // set-up query for derived variable A
        InterpolationQueryParticle query_derived(n_local);
        query_derived.setCoords(0, interp_x_local.data())
            .setCoords(1, interp_y_local.data())
            .setCoords(2, interp_z_local.data())
            .addComp(0, A_local.data(), VariableType::derived, BCParity::even,
                     Derivative::LOCAL);

        // set-up query for state variable B
        InterpolationQueryParticle query_state(n_local);
        query_state.setCoords(0, interp_x_local.data())
            .setCoords(1, interp_y_local.data())
            .setCoords(2, interp_z_local.data())
            .addComp(c_polystate, B_local.data(), VariableType::state);

        // set up interpolation using Particles for derived vars
        ParticleInterpolator<1> interpolator_derived;

        interpolator_derived.setup(&gr_amr);
        interpolator_derived.interp(
            query_derived, false, PolynomialDerivedQuantity::name,
            0.0); // do not refresh particles as the query remains the same

        // set up interpolation using Particles for state vars
        ParticleInterpolator<1> interpolator_state;
        interpolator_state.setup(&gr_amr);
        interpolator_state.interp(
            query_state,
            false); // do not refresh particles as the query remains the same

        for (int ipoint = 0; ipoint < n_local; ++ipoint)
        {
            double x = interp_x_local[ipoint] - center[0];
            double y = interp_y_local[ipoint] - center[1];
            double z = interp_z_local[ipoint] - center[2];

            double A_known = 42. + x * x + y * y * z * z;
            double B_known = pow(z, 3);

            INFO("Interpolated A is "
                 << A_local[ipoint] << " at point x = " << x << " y = " << y
                 << " z = " << z << ". The true value should be " << A_known);
            INFO("Interpolated B is "
                 << B_local[ipoint] << " at point x = " << x << " y = " << y
                 << " z = " << z << ". The true value should be " << B_known);

            CHECK(A_local[ipoint] == doctest::Approx(A_known).epsilon(1e-10));
            CHECK(B_local[ipoint] == doctest::Approx(B_known).epsilon(1e-10));
        }
    }
    amrex::Finalize();
}
