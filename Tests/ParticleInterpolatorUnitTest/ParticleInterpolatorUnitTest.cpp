/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// Doctest header
#include "doctest.h"

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
        this_file.parent_path() /
        std::filesystem::path("ParticleInterpolatorUnitTest.inputs");
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
        ParticleInterpolatorLevel::variableSetUp();

        // Set the center
        PolynomialDerivedQuantity::set_center(sim_params.center);

        // Set up the AMR object
        DefaultLevelFactory<ParticleInterpolatorLevel>
            interpolator_test_level_fact;
        GRAMR gr_amr(&interpolator_test_level_fact);
        gr_amr.init(0., sim_params.stop_time);

        // Build the point from sim_params
        const int num_points = sim_params.num_points;

        std::vector<double> A(num_points); // for storing derived polynomial
        std::vector<double> B(num_points); // for storing state polynomial
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

        // set-up query for derived variable A
        InterpolationQueryParticle query(num_points);
        query.setCoords(0, interp_x.data())
            .setCoords(1, interp_y.data())
            .setCoords(2, interp_z.data())
            .addComp(0, A.data(), Derivative::LOCAL, VariableType::derived);

        // set up interpolation using Particles for derived vars
        ParticleInterpolator<1> interpolator;
        interpolator.set_gramr_ptr(&gr_amr, sim_params.boundary_params, 0,
                                   true);
        interpolator.set_derived_var_parity(0, BCParity::even);
        int ngrow = 2;
        interpolator.interp(query, VariableType::derived,
                            PolynomialDerivedQuantity::name, 0.0);

        // set-up query for state variable B
        InterpolationQueryParticle query_state(num_points);
        query_state.setCoords(0, interp_x.data())
            .setCoords(1, interp_y.data())
            .setCoords(2, interp_z.data())
            .addComp(0, B.data(), Derivative::LOCAL, VariableType::state);

        // set up interpolation using Particles for state vars
        ParticleInterpolator<1> interpolator_state;
        interpolator_state.set_gramr_ptr(&gr_amr, sim_params.boundary_params, 0,
                                         true);
        interpolator_state.interp(query_state, VariableType::state);

        if (amrex::ParallelDescriptor::MyProc() == 0)
        {
            for (int ipoint = 0; ipoint < num_points; ++ipoint)
            {
                double x = interp_x[ipoint] - sim_params.center[0];
                double y = interp_y[ipoint] - sim_params.center[1];
                double z = interp_z[ipoint] - sim_params.center[2];

                double A_known = 42. + x * x + y * y * z * z; // derived
                double B_known = pow(x, 3);                   // state

                INFO("Interpolated A is " << A[ipoint] << " at point x = " << x
                                          << " y = " << y << " z = " << z
                                          << ". The true value should be "
                                          << A_known);
                INFO("Interpolated B is " << B[ipoint] << " at point x = " << x
                                          << " y = " << y << " z = " << z
                                          << ". The true value should be "
                                          << B_known);

                CHECK(A[ipoint] == doctest::Approx(A_known).epsilon(1e-10));
                CHECK(B[ipoint] == doctest::Approx(B_known).epsilon(1e-10));
            }
        }
    }

    amrex::Finalize();
}
