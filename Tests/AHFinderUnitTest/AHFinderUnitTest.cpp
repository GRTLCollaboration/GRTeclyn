// Test include
#include "AHFinderUnitTest.hpp"

// Doctest header
#include "doctest.h"

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
#include "DefaultLevelBld.hpp"
#include "GRAmr.hpp"
#include "GRAmrLevel.hpp"
#include "GRParmParse.hpp"

// Problem specific includes
#include "AHFinder.hpp"
#include "AHFinderLevel.hpp"
#include "BoostedBHInitialData.hpp"
#include "SimulationParameters.hpp"

// Others
#include <filesystem>

void run_ah_finder_unit_test()
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
        AHFinderLevel::variableSetUp();

        // Set up the AMR object
        DefaultLevelBld<AHFinderLevel> ahfinder_test_level_fact;
        GRAmr gr_amr(&ahfinder_test_level_fact);

        amrex::Real stop_time{};
        pp.get("evolution.stop_time", stop_time);
        gr_amr.init(0., stop_time);

        // Read from params
        int num_particles{};
        pp.get("test.num_particles", num_particles);

        bool verbosity{};
        pp.get("particle_interpolator.verbosity", verbosity);

        // Number of processes and local processes
        const int nprocs = amrex::ParallelDescriptor::NProcs();
        const int myproc = amrex::ParallelDescriptor::MyProc();

        // Partition the points across ranks
        const int base      = num_particles / nprocs;
        const int remainder = num_particles % nprocs;

        const int n_local =
            base + (myproc < remainder
                        ? 1
                        : 0); // local number of particles on the rank
        if (verbosity)
        {
            amrex::AllPrint() << "I am rank " << myproc << " and I query "
                              << n_local << " LOCAL particles \n";
        }

        // By default, search for the individual horizon around puncture A.
        // test.ah_offset / test.ah_guess_radius override this, which is how the
        // common horizon of a close binary is targeted instead -- see
        // docs/ah_finder_implicit_ptc.md.
        BoostedBHInitialData::params_t bh1_params(1);
        bh1_params.fill_params();

        std::array<double, AMREX_SPACEDIM> ah_center = bh1_params.center;
        amrex::Real guess_radius                     = 0.5 * bh1_params.mass;

        // test.ah_offset is measured from the domain centre, the same
        // convention as bh1.offset / bh2.offset, so (0,0,0) is the midpoint of
        // the binary and targets the common horizon.
        std::vector<amrex::Real> ah_offset;
        if (pp.queryarr("test.ah_offset", ah_offset) &&
            ah_offset.size() == AMREX_SPACEDIM)
        {
            GRParmParse geom_pp("geometry");
            std::array<amrex::Real, AMREX_SPACEDIM> domain_center{};
            geom_pp.get("center", domain_center);

            for (int i = 0; i < AMREX_SPACEDIM; ++i)
            {
                ah_center[i] = domain_center[i] + ah_offset[i];
            }
        }
        pp.query("test.ah_guess_radius", guess_radius);

        AHFinder<21> finder(num_particles, ah_center, guess_radius);

        finder.init(&gr_amr);
        finder.find();
    }

    amrex::Finalize();
}
