/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SETUP_FUNCTIONS_HPP_
#define SETUP_FUNCTIONS_HPP_
// This file incldues several functions that need to be called to
// set up the runs but aren't very interesting for the normal user.

// Other includes
#include <iostream>
using std::cerr;
using std::endl;
#include "AMReXParameters.hpp"
#include "DerivativeSetup.hpp"
#include "FilesystemTools.hpp"
#include "GRAMR.hpp"
#include "GRParmParse.hpp"
#include "IntegrationMethodSetup.hpp"

#include "simd.hpp"

#ifdef EQUATION_DEBUG_MODE
#include "DebuggingTools.hpp"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

/// This function calls MPI_Init, makes sure a parameter file is supplied etc...
void mainSetup(int argc, char *argv[]);

/// This function calls all finalisations
void mainFinalize();

/// Sets up the grid parameters, problem domain and AMR object
//xxxxxvoid setupAMRObject(AMR &gr_amr, AMRLevelFactory &a_factory);

void mainSetup(int argc, char *argv[])
{
#ifdef CH_MPI
    // Start MPI
    MPI_Init(&argc, &argv);
#ifdef CH_AIX
    H5dont_atexit();
#endif
// setChomboMPIErrorHandler();
#endif

    int rank, number_procs;
#ifdef CH_MPI
    MPI_Comm_rank(Chombo_MPI::comm, &rank);
    MPI_Comm_size(Chombo_MPI::comm, &number_procs);
#else
    rank = 0;
    number_procs = 1;
#endif

#ifdef EQUATION_DEBUG_MODE
    EquationDebugging::check_no_omp();
    amrex::Warning("GRChombo is running in equation debug mode. This mode is "
                    "intended only for debugging and leads to significantly "
                    "worse performance.");
#endif

    if (rank == 0)
    {
        std::cout << " number_procs = " << number_procs << endl;
#ifdef _OPENMP
        std::cout << " threads = " << omp_get_max_threads() << endl;
#endif
        std::cout << " simd width (doubles) = " << simd_traits<double>::simd_len
                  << endl;
    }

    const int required_argc = 2;
    if (argc < required_argc)
    {
        cerr << " usage " << argv[0] << " <input_file_name> " << endl;
        exit(0);
    }
}

void mainFinalize()
{
#ifdef CH_MPI
    // Exit MPI
#ifdef CH_USE_MEMORY_TRACKING
    dumpmemoryatexit();
#endif
    MPI_Finalize();
#endif
}

#if 0
//xxxxx
void setupAMRObject(GRAMR &gr_amr, AMRLevelFactory &a_factory)
{
    // Reread the params - just the base ones
    // Note that we could have passed these through the function
    // but this way preserves backwards compatibility
    GRParmParse pp;
    AMReXParameters amrex_params(pp);

    // set size of box
    Box problem_domain(IntVect::Zero, amrex_params.ivN);
    ProblemDomain physdomain(problem_domain);

    // set periodicity
    for (int dir = 0; dir < SpaceDim; dir++)
    {
        physdomain.setPeriodic(dir,
                               amrex_params.boundary_params.is_periodic[dir]);
    }

    // Define the AMR object
    gr_amr.define(amrex_params.max_level, amrex_params.ref_ratios, physdomain,
                  &a_factory);

    // The buffer defines the minimum number of level l cells there have to be
    // between level l+1 and level l-1
    // It needs to be at least ceil(num_ghosts/max_ref_ratio) for proper nesting
    gr_amr.gridBufferSize(amrex_params.grid_buffer_size);

    // set checkpoint and plot intervals and prefixes
#ifdef CH_USE_HDF5
    gr_amr.checkpointInterval(amrex_params.checkpoint_interval);
    gr_amr.checkpointPrefix(amrex_params.hdf5_path +
                            amrex_params.checkpoint_prefix);
    if (amrex_params.plot_interval != 0)
    {
        gr_amr.plotInterval(amrex_params.plot_interval);
        gr_amr.plotPrefix(amrex_params.hdf5_path + amrex_params.plot_prefix);
    }
#endif

    // Number of coarse time steps from one regridding to the next
    gr_amr.regridIntervals(amrex_params.regrid_interval);

    // max and min box sizes, fill ratio determining accuracy of regrid
    gr_amr.maxGridSize(amrex_params.max_grid_size);
    gr_amr.blockFactor(amrex_params.block_factor);
    gr_amr.fillRatio(amrex_params.fill_ratio);

    // Set verbosity
    gr_amr.verbosity(amrex_params.verbosity);

    // Set timeEps to half of finest level dt
    // Chombo sets it to 1.e-6 by default (AMR::setDefaultValues in AMR.cpp)
    // This is only not enough for >~20 levels
    double eps = 1.;
    for (int ilevel = 0; ilevel < amrex_params.max_level; ++ilevel)
        eps /= amrex_params.ref_ratios[ilevel];
    gr_amr.timeEps(std::min(1.e-6, eps / 2.));

    // Set up input files
    if (!amrex_params.restart_from_checkpoint)
    {
#ifdef CH_USE_HDF5
        if (!FilesystemTools::directory_exists(amrex_params.hdf5_path))
            FilesystemTools::mkdir_recursive(amrex_params.hdf5_path);
#endif

        gr_amr.setupForNewAMRRun();
    }
    else
    {
#ifdef CH_USE_HDF5
        HDF5Handle handle(amrex_params.restart_file, HDF5Handle::OPEN_RDONLY);
        // read from checkpoint file
        gr_amr.setupForRestart(handle);
        handle.close();
#else
        amrex::Abort("GRChombo restart only defined with hdf5");
#endif
    }
}
#endif

#endif /* SETUP_FUNCTIONS_HPP_ */
