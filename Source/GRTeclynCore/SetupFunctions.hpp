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
#include "GRAmr.hpp"
#include "GRParmParse.hpp"
#include "IntegrationMethodSetup.hpp"
#include "SimulationParameters.hpp"
#ifdef USE_STOIC_QUOTES
#include "StoicQuotes.hpp"
#endif

#ifdef EQUATION_DEBUG_MODE
#include "DebuggingTools.hpp"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include <fstream>
#include <iostream>
#ifdef USE_STOIC_QUOTES
#include <string>
#endif

#ifndef GRTECLYN_VERSION
#define GRTECLYN_VERSION "unknown"
#endif

/// Return true when the GRTeclyn-only flag following the AMReX `--` separator
/// requests parameter validation without running the simulation.
[[nodiscard]] inline bool just_check_params()
{
    // AMReX reserves argument 0 for the executable name.
    for (int argument_index = 1;
         argument_index <= amrex::command_argument_count(); ++argument_index)
    {
        if (amrex::get_command_argument(argument_index) == "-check_params")
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

/// Print the job status followed by a randomly selected Stoic quote.
void print_job_end_message(int status);

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
        GRParmParse::prettyPrintTable(ofs);
    }

    amrex::Print() << "GRTeclyn (" << GRTECLYN_VERSION << ") initialized"
                   << '\n';

    if (just_check_params())
    {
        if (GRParmParse::warnings_issued())
        {
            amrex::Print()
                << "User-specified parameter checks were applied and some "
                   "warnings were identified, which should be listed above. "
                   "You can use parameters_and_version.txt to check that the "
                   "parameters were set as you intended. Warnings won't "
                   "prevent the code from running, but you should check that "
                   "you understand why they can be safely ignored before "
                   "continuing.\n";
        }
        else
        {
            amrex::Print()
                << "User-specified parameter checks were applied and no "
                   "warnings or errors were identified. You can use "
                   "parameters_and_version.txt to check that the parameters "
                   "were set as you intended. This does not guarantee that "
                   "your simulation will work, but it is a good start.\n";
        }
    }

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

void print_job_end_message(int status)
{
    if (status == 0)
    {
        amrex::Print() << "GRTeclyn finished.\n";
    }
    else
    {
        amrex::Print() << "GRTeclyn failed with return code " << status << '\n';
    }

#ifdef USE_STOIC_QUOTES
    {
        const std::string message =
            " Stoic wisdom: " +
            std::string(StoicQuotes::random_quote(status == 0)) + " ";
        const std::string border(message.size(), '-');
        amrex::Print() << '+' << border << "+\n"
                       << '|' << message << "|\n"
                       << '+' << border << "+\n";
    }
#endif
}

#endif /* SETUP_FUNCTIONS_HPP_ */
