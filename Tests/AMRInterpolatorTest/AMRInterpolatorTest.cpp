// Catch2 header
#include "catch_amalgamated.hpp"

// Common test headers
#include "TestsArgs.hpp"

// General includes:
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>

// GRTeclyn headers
#include "DefaultLevelFactory.hpp"
#include "GRAMR.hpp"
#include "GRParmParse.hpp"
#include "SetupFunctions.hpp"
#include "SimulationParameters.hpp"

// Problem specific includes:
// #include "AMRInterpolator.hpp"
// #include "InterpolationQuery.hpp"
#include "InterpolatorTestLevel.hpp"
// #include "Lagrange.hpp"
#include "UserVariables.hpp"

TEST_CASE("AMRInterpolator")
{
    // Fabricate argc and argv
    int fake_argc       = 2;
    char *exename       = Tests::g_args.get_argv()[0];
    char param_file[50] = "./AMRInterpolatorTest/AMRInterpolatorTest.inputs";

    char *fake_argv[fake_argc];
    fake_argv[0] = exename;
    fake_argv[1] = param_file;

    mainSetup(fake_argc, fake_argv);
    {
        GRParmParse pp;
        SimulationParameters sim_params(pp);

        GRAMR::set_simulation_parameters(sim_params);

        DefaultLevelFactory<InterpolatorTestLevel> interpolator_test_level_fact;

        GRAMR gr_amr(&interpolator_test_level_fact);

        // Setup the AMRInterpolator
        const int num_points = sim_params.num_points;

        std::vector<double> A(num_points);
        std::vector<double> B(num_points);
        std::vector<double> B_dx(num_points);
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

        // InterpolationQuery query(num_points);
        // query.setCoords(0, interp_x.data())
        //     .setCoords(1, interp_y.data())
        //     .setCoords(2, interp_z.data())
        //     .addComp(c_A, A.data())
        //     .addComp(c_B, B.data())
        // .addComp(c_B, B_dx.data(), Derivative::dx);

        // AMRInterpolator<Lagrange<4>> interpolator(gr_amr, sim_params.origin,
        //                                           sim_params.dx,
        //                                           sim_params.boundary_params,
        //                                           0);
        // interpolator.interp(query);

        for (int ipoint = 0; ipoint < num_points; ++ipoint)
        {
            double x = interp_x[ipoint] - sim_params.center[0];
            double y = interp_y[ipoint] - sim_params.center[1];
            double z = interp_z[ipoint] - sim_params.center[2];

            double value_A    = 42. + x * x + y * y * z * z;
            double value_B    = std::pow(x, 3);
            double value_B_dx = 3. * std::pow(x, 2);

            CHECK_THAT(A[ipoint], Catch::Matchers::WithinAbs(value_A, 1e-10));
            CHECK_THAT(B[ipoint], Catch::Matchers::WithinAbs(value_B, 1e-10));
            CHECK_THAT(B_dx[ipoint],
                       Catch::Matchers::WithinAbs(value_B_dx, 1e-10));
        }
    }
    mainFinalize();
}