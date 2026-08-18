/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SETUP_FUNCTIONS_HPP_
#define SETUP_FUNCTIONS_HPP_
// This file incldues several functions that need to be called to
// set up the runs but aren't very interesting for the normal user.

// xxxxx various setups
#include "DerivativeSetup.hpp"
#include "FilesystemTools.hpp"
#include "GRAMR.hpp"
#include "GRParmParse.hpp"
#include "IntegrationMethodSetup.hpp"
#include "SimulationParameters.hpp"

#ifdef EQUATION_DEBUG_MODE
#include "DebuggingTools.hpp"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include <fstream>
#include <iostream>
#include <string_view>

#ifndef GRTECLYN_VERSION
#define GRTECLYN_VERSION "unknown"
#endif

// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
inline bool just_check_params(int argc, char *argv[])
{
    for (int i = 1; i < argc; ++i)
    {
        if (std::string_view(argv[i]) == "-check_params")
        {
            return true;
        }
    }
    return false;
}

// NOLINTBEGIN(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
/// This function calls MPI_Init
void mainSetup(int argc, char *argv[]);

/// This function calls all finalisations
void mainFinalize();

void mainSetup(int argc, char *argv[])
{
    // NOLINTEND(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
    // NOLINTNEXTLINE(bugprone-casting-through-void)
    amrex::Initialize(
        argc, argv, std::function<void()>(SimulationParameters::check_params));

    if (amrex::ParallelDescriptor::IOProcessor())
    {
        std::ofstream ofs("parameters_and_version.txt");
        ofs << "Using GRTeclyn version (" << GRTECLYN_VERSION << ")\n\n";
        GRParmParse pp;
        pp.prettyPrintTable(ofs);
    }

    amrex::Print() << "GRTeclyn (" << GRTECLYN_VERSION << ") initialized"
                   << '\n';

#ifdef EQUATION_DEBUG_MODE
    EquationDebugging::check_no_omp();
    amrex::Warning("GRTeclyn is running in equation debug mode. This mode is "
                   "intended only for debugging and leads to significantly "
                   "worse performance.");
#endif

    const int required_argc = 2;
    if (argc < required_argc)
    {
        amrex::Finalize();
        std::cerr << " usage " << argv[0] << " <input_file_name> " << '\n';
        exit(0);
    }
}

void mainFinalize() { amrex::Finalize(); }

#endif /* SETUP_FUNCTIONS_HPP_ */
