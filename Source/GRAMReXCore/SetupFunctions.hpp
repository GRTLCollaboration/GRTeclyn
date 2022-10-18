/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SETUP_FUNCTIONS_HPP_
#define SETUP_FUNCTIONS_HPP_
// This file incldues several functions that need to be called to
// set up the runs but aren't very interesting for the normal user.

//xxxxx various setups
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

#include <iostream>

/// This function calls MPI_Init, makes sure a parameter file is supplied etc...
void mainSetup(int argc, char *argv[]);

/// This function calls all finalisations
void mainFinalize();

#if !defined(AMREX_USE_GPU)
const int simd_traits<double>::simd_len; // Still needs to be defined
#endif

void mainSetup(int argc, char *argv[])
{
    bool use_parm_parse = true;
    amrex::Initialize(argc, argv, use_parm_parse, MPI_COMM_WORLD,
                      [] () {
                          amrex::ParmParse pp("amrex");
                          // don't use managed memory
                          pp.add("the_arena_is_managed", false);
                      });

#ifdef EQUATION_DEBUG_MODE
    EquationDebugging::check_no_omp();
    amrex::Warning("GRAMReX is running in equation debug mode. This mode is "
                   "intended only for debugging and leads to significantly "
                   "worse performance.");
#endif

#if !defined(AMREX_USE_GPU)
    amrex::Print() << " simd width (doubles) = " << simd_traits<double>::simd_len
                   << std::endl;
#endif

    const int required_argc = 2;
    if (argc < required_argc)
    {
        amrex::Finalize();
        std::cerr << " usage " << argv[0] << " <input_file_name> " << std::endl;
        exit(0);
    }
}

void mainFinalize()
{
    amrex::Finalize();
}

#endif /* SETUP_FUNCTIONS_HPP_ */
