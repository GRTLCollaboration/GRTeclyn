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
#include "DefaultLevelFactory.hpp"
#include "GRAMR.hpp"
#include "GRAMRLevel.hpp"
#include "GRParmParse.hpp"

// Problem specific includes
#include "AHFinder.hpp"
#include "AHFinderLevel.hpp"
#include "SimulationParameters.hpp"

// Others
#include <filesystem>

void run_ah_finder_unit_test()
{
    // Use an input file that is in the same directory as this file for the
    // second argument
    std::filesystem::path this_file(__FILE__);
    std::filesystem::path input_file =
        this_file.parent_path() /
        std::filesystem::path("AHFinderUnitTest.inputs");
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
        AHFinderLevel::variableSetUp();

        // Set up the AMR object
        DefaultLevelFactory<AHFinderLevel> ahfinder_test_level_fact;
        GRAMR gr_amr(&ahfinder_test_level_fact);
        gr_amr.init(0., sim_params.stop_time);

        // Read from params
        const int num_particles = 64;
        bool verbosity          = sim_params.verbosity;

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

        AHFinder<1> finder(num_particles, sim_params.center);

        finder.init(&gr_amr, sim_params.boundary_params, verbosity);
    }

    amrex::Finalize();
}