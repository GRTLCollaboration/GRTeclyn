/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// Doctest header
#include "doctest.h"

// State vars
#include "StateVariables.hpp"

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
    amrex::Initialize(new_argc, new_argv);
    {
        // Load the parameter file and construct SimulationParameters
        GRParmParse pp;
        SimulationParameters sim_params(pp);
        GRAMR::set_simulation_parameters(sim_params);
        SphericalExtractionTestLevel::variableSetUp();

        DefaultLevelFactory<SphericalExtractionTestLevel>
            surface_extraction_test_level_fact;
        GRAMR gr_amr(&surface_extraction_test_level_fact);

        gr_amr.init(0., sim_params.stop_time);

        bool broadcast_integral = true;

        std::pair<std::vector<amrex::ParticleReal>,
                  std::vector<amrex::ParticleReal>>
            integral_lo_trapezium, integral_hi_trapezium;
        std::pair<std::vector<amrex::ParticleReal>,
                  std::vector<amrex::ParticleReal>>
            integral_lo_simpson, integral_hi_simpson;
        std::pair<std::vector<amrex::ParticleReal>,
                  std::vector<amrex::ParticleReal>>
            integral_lo_boole, integral_hi_boole;

        std::vector<int> state_vars = {c_phi_Re, c_phi_Im};

        // Real part is the zeroth component and imaginary part is the first
        // component
        SphericalExtraction<2>::complex_function_t extracted_harmonic =
            [](std::vector<amrex::ParticleReal> &data, amrex::ParticleReal,
               amrex::ParticleReal, amrex::ParticleReal)
        { return std::make_pair(data[0], data[1]); };

        {
            // Initiate ParticleInterpolator
            ParticleInterpolator<2> interpolator;
            interpolator.setup(&gr_amr, sim_params.boundary_params,
                               sim_params.verbosity);
            // Low resolution spherical extraction
            SphericalExtraction<2> spherical_extraction_lo(
                sim_params.extraction_params_lo, state_vars,
                sim_params.coarsest_dx * sim_params.dt_multiplier, 0.0, true,
                0.0);

            spherical_extraction_lo.extract(&interpolator);
            spherical_extraction_lo.write_extraction("ExtractionOutLo_");

            // Add the spherical harmonic mode integrands for each resolution
            // and for the trapezium rule, Simpson's rule and Boole's rule
            // Always use trapezium rule in phi as this is periodic
            spherical_extraction_lo.add_mode_integrand(
                sim_params.es, sim_params.el, sim_params.em, extracted_harmonic,
                integral_lo_trapezium, IntegrationMethod::trapezium,
                IntegrationMethod::trapezium, broadcast_integral);
            spherical_extraction_lo.add_mode_integrand(
                sim_params.es, sim_params.el, sim_params.em, extracted_harmonic,
                integral_lo_simpson, IntegrationMethod::simpson,
                IntegrationMethod::trapezium, broadcast_integral);
            spherical_extraction_lo.add_mode_integrand(
                sim_params.es, sim_params.el, sim_params.em, extracted_harmonic,
                integral_lo_boole, IntegrationMethod::boole,
                IntegrationMethod::trapezium, broadcast_integral);

            // Do the surface integration
            spherical_extraction_lo.integrate();
        }

        {
            ParticleInterpolator<2> interpolator_hi;
            interpolator_hi.setup(&gr_amr, sim_params.boundary_params,
                                  sim_params.verbosity);

            // High resolution spherical extraction
            spherical_extraction_params_t extraction_params_hi =
                sim_params.extraction_params_lo;

            // We are only checking the convergence in theta integration
            extraction_params_hi.num_points_theta() *= 2;
            // Need to subtract a point as it's the number of subintervals we
            // want to double for theta
            extraction_params_hi.num_points_theta() -= 1;

            SphericalExtraction<2> spherical_extraction_hi(
                extraction_params_hi, state_vars,
                sim_params.coarsest_dx * sim_params.dt_multiplier, 0.0, true,
                0.0);

            spherical_extraction_hi.extract(&interpolator_hi);
            spherical_extraction_hi.write_extraction("ExtractionOutHi_");
            spherical_extraction_hi.add_mode_integrand(
                sim_params.es, sim_params.el, sim_params.em, extracted_harmonic,
                integral_hi_trapezium, IntegrationMethod::trapezium,
                IntegrationMethod::trapezium, broadcast_integral);
            spherical_extraction_hi.add_mode_integrand(
                sim_params.es, sim_params.el, sim_params.em, extracted_harmonic,
                integral_hi_simpson, IntegrationMethod::simpson,
                IntegrationMethod::trapezium, broadcast_integral);
            spherical_extraction_hi.add_mode_integrand(
                sim_params.es, sim_params.el, sim_params.em, extracted_harmonic,
                integral_hi_boole, IntegrationMethod::boole,
                IntegrationMethod::trapezium, broadcast_integral);
            spherical_extraction_hi.integrate();
        }

        amrex::Print() << std::setprecision(10);

        for (int iradius = 0;
             iradius < sim_params.extraction_params_lo.num_extraction_radii();
             ++iradius)
        {
            amrex::ParticleReal r =
                sim_params.extraction_params_lo.extraction_radii()[iradius];

            // NOLINTBEGIN(cppcoreguidelines-init-variables)
            amrex::ParticleReal integral_re_lo_trapezium =
                integral_lo_trapezium.first[iradius];
            amrex::ParticleReal integral_re_hi_trapezium =
                integral_hi_trapezium.first[iradius];
            amrex::ParticleReal integral_re_lo_simpson =
                integral_lo_simpson.first[iradius];
            amrex::ParticleReal integral_re_hi_simpson =
                integral_hi_simpson.first[iradius];
            amrex::ParticleReal integral_re_lo_boole =
                integral_lo_boole.first[iradius];
            amrex::ParticleReal integral_re_hi_boole =
                integral_hi_boole.first[iradius];
            // NOLINTEND(cppcoreguidelines-init-variables)

            amrex::ParticleReal analytic_integral = 1.0;

            amrex::ParticleReal convergence_factor_trapezium =
                std::abs((integral_re_lo_trapezium - analytic_integral) /
                         (integral_re_hi_trapezium - analytic_integral));
            amrex::ParticleReal convergence_factor_simpson =
                std::abs((integral_re_lo_simpson - analytic_integral) /
                         (integral_re_hi_simpson - analytic_integral));
            amrex::ParticleReal convergence_factor_boole =
                std::abs((integral_re_lo_boole - analytic_integral) /
                         (integral_re_hi_boole - analytic_integral));

            amrex::ParticleReal convergence_order_trapezium =
                std::log2(convergence_factor_trapezium);
            amrex::ParticleReal convergence_order_simpson =
                std::log2(convergence_factor_simpson);
            amrex::ParticleReal convergence_order_boole =
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