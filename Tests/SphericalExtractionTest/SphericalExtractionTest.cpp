/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// Doctest header
#include "doctest.h"

// Common includes
#include "doctestCLIArgs.hpp"

// AMReX includes
#include "AMReX.H"
#include <AMReX_Print.H>

// Base includes
#include "DefaultLevelFactory.hpp"
#include "GRAMR.hpp"
#include "GRParmParse.hpp"

// Problem specific includes
#include "IntegrationMethod.hpp"
#include "ParticleInterpolator.hpp"
#include "SimulationParameters.hpp"
#include "SphericalExtraction.hpp"
#include "SphericalExtractionTestLevel.hpp"

// Others
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <utility>
#include <vector>

void run_spherical_extraction_test()
{
    // Use an input file that is in the same directory as this file
    std::filesystem::path this_file(__FILE__);
    std::filesystem::path input_file =
        this_file.parent_path() /
        std::filesystem::path("SphericalExtractionTest.inputs");
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
        // Load the parameter file and construct SimulationParameters
        GRParmParse pp;
        SphericalExtractionTestLevel::variableSetUp();

        DefaultLevelFactory<SphericalExtractionTestLevel>
            surface_extraction_test_level_fact;
        GRAMR gr_amr(&surface_extraction_test_level_fact);

        double stop_time{};
        pp.get("grteclyn.stop_time", stop_time);
        gr_amr.init(0., stop_time);

        double coarsest_dx;
        double dt_multiplier;

        pp.get("grteclyn.dt_multiplier", dt_multiplier);
        pp.get("grteclyn.coarsest_dx", coarsest_dx);
        bool broadcast_integral = true;

        std::pair<std::vector<double>, std::vector<double>>
            integral_lo_trapezium, integral_hi_trapezium;
        std::pair<std::vector<double>, std::vector<double>> integral_lo_simpson,
            integral_hi_simpson;
        std::pair<std::vector<double>, std::vector<double>> integral_lo_boole,
            integral_hi_boole;

        std::vector<int> state_vars = {c_phi, c_Pi};

        // Real part is the zeroth component and imaginary part is the first
        // component
        SphericalExtraction<2>::complex_function_t extracted_harmonic =
            [](std::vector<double> &data, double, double, double)
        { return std::make_pair(data[0], data[1]); };
        int es{0};
        int el{2};
        int em{0}; // spherical harmonic params
        pp.queryAdd("extraction.es", es);
        pp.queryAdd("extraction.el", el);
        pp.queryAdd("extraction.em", em);

        int num_points_theta_lo{};
        pp.get("extraction.num_points_theta_lo", num_points_theta_lo);
        pp.add("extraction.num_points_theta", num_points_theta_lo);
        {
            // Initiate ParticleInterpolator
            ParticleInterpolator<2> interpolator;
            interpolator.setup(&gr_amr);
            // Low resolution spherical extraction
            spherical_extraction_params_t extraction_params;
            extraction_params.fill_params();
            SphericalExtraction<2> spherical_extraction_lo(
                extraction_params, state_vars, coarsest_dx * dt_multiplier, 0.0,
                true, 0.0);

            spherical_extraction_lo.extract(&interpolator);
            spherical_extraction_lo.write_extraction("ExtractionOutLo_");

            // Add the spherical harmonic mode integrands for each resolution
            // and for the trapezium rule, Simpson's rule and Boole's rule
            // Always use trapezium rule in phi as this is periodic
            spherical_extraction_lo.add_mode_integrand(
                es, el, em, extracted_harmonic, integral_lo_trapezium,
                IntegrationMethod::trapezium, IntegrationMethod::trapezium,
                broadcast_integral);
            spherical_extraction_lo.add_mode_integrand(
                es, el, em, extracted_harmonic, integral_lo_simpson,
                IntegrationMethod::simpson, IntegrationMethod::trapezium,
                broadcast_integral);
            spherical_extraction_lo.add_mode_integrand(
                es, el, em, extracted_harmonic, integral_lo_boole,
                IntegrationMethod::boole, IntegrationMethod::trapezium,
                broadcast_integral);

            // Do the surface integration
            spherical_extraction_lo.integrate();
        }

        {
            ParticleInterpolator<2> interpolator_hi;
            interpolator_hi.setup(&gr_amr);

            // We are only checking the convergence in theta integration

            int num_points_theta_hi  = num_points_theta_lo;
            num_points_theta_hi     *= 2;
            // Need to subtract a point as it's the number of subintervals we
            // want to double for theta
            num_points_theta_hi -= 1;
            pp.add("extraction.num_points_theta", num_points_theta_hi);

            spherical_extraction_params_t extraction_params;
            extraction_params.fill_params();
            SphericalExtraction<2> spherical_extraction_hi(
                extraction_params, state_vars, coarsest_dx * dt_multiplier, 0.0,
                true, 0.0);

            spherical_extraction_hi.extract(&interpolator_hi);
            spherical_extraction_hi.write_extraction("ExtractionOutHi_");
            spherical_extraction_hi.add_mode_integrand(
                es, el, em, extracted_harmonic, integral_hi_trapezium,
                IntegrationMethod::trapezium, IntegrationMethod::trapezium,
                broadcast_integral);
            spherical_extraction_hi.add_mode_integrand(
                es, el, em, extracted_harmonic, integral_hi_simpson,
                IntegrationMethod::simpson, IntegrationMethod::trapezium,
                broadcast_integral);
            spherical_extraction_hi.add_mode_integrand(
                es, el, em, extracted_harmonic, integral_hi_boole,
                IntegrationMethod::boole, IntegrationMethod::trapezium,
                broadcast_integral);
            spherical_extraction_hi.integrate();
        }

        amrex::Print() << std::setprecision(10);

        int num_extraction_radii{};
        pp.get("extraction.num_extraction_radii", num_extraction_radii);
        std::vector<int> extraction_radii_stdvect(num_extraction_radii);
        pp.queryAdd("extraction.extraction_radii", extraction_radii_stdvect);

        for (int iradius = 0; iradius < num_extraction_radii; ++iradius)
        {
            double r = extraction_radii_stdvect[iradius];

            // NOLINTBEGIN(cppcoreguidelines-init-variables)
            double integral_re_lo_trapezium =
                integral_lo_trapezium.first[iradius];
            double integral_re_hi_trapezium =
                integral_hi_trapezium.first[iradius];
            double integral_re_lo_simpson = integral_lo_simpson.first[iradius];
            double integral_re_hi_simpson = integral_hi_simpson.first[iradius];
            double integral_re_lo_boole   = integral_lo_boole.first[iradius];
            double integral_re_hi_boole   = integral_hi_boole.first[iradius];
            // NOLINTEND(cppcoreguidelines-init-variables)

            double analytic_integral = 1.0;

            double convergence_factor_trapezium =
                std::abs((integral_re_lo_trapezium - analytic_integral) /
                         (integral_re_hi_trapezium - analytic_integral));
            double convergence_factor_simpson =
                std::abs((integral_re_lo_simpson - analytic_integral) /
                         (integral_re_hi_simpson - analytic_integral));
            double convergence_factor_boole =
                std::abs((integral_re_lo_boole - analytic_integral) /
                         (integral_re_hi_boole - analytic_integral));

            double convergence_order_trapezium =
                std::log2(convergence_factor_trapezium);
            double convergence_order_simpson =
                std::log2(convergence_factor_simpson);
            double convergence_order_boole =
                std::log2(convergence_factor_boole);

            INFO("At r = " << r);

            if (amrex::ParallelDescriptor::MyProc() == 0)
            {
                // Trapezium rule should have second order convergence
                CHECK_GT(convergence_order_trapezium, 1.5);

                // Simpson's rule should have fourth order convergence
                CHECK_GT(convergence_order_simpson, 3.5);

                // Boole's rule should have sixth order convergence
                CHECK_GT(convergence_order_boole, 5.5);
            }
        }
    }
    amrex::Finalize();
}