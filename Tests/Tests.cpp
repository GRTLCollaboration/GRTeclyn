// Doctest header
#ifdef AMREX_USE_SYCL
// Intel's GPU runtime uses SIGSEGV to trigger migration of managed memory
#define DOCTEST_CONFIG_NO_POSIX_SIGNALS
#endif
#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"

#include "TestCases.hpp" // Test cases are defined here
#include "doctestOutput.hpp"

// system headers
#include <iomanip>
#include <iostream>

#include "AMReX.H"
#include "AMReX_REAL.H"
#include "AMReX_ccse-mpi.H"

int main(int argc, char *argv[])
{
#ifdef BL_USE_MPI
    // We can only initialize and finalize MPI once so do it here
    MPI_Init(&argc, &argv);
#endif

    doctest::Context doctest_context(argc, argv);
#ifdef BL_USE_MPI
    doctest_context.setCout(
        &doctest::hide_output_from_non_zero_ranks(std::cout));
#endif

    // Default AMReX verbosity to 0 to avoid the "Initialized", "Finalized" and
    // memory usage messages
    amrex::system::verbose = 0;

    constexpr int cout_precision = 17;
    std::cout << std::setprecision(cout_precision);

    int result = doctest_context.run();

#ifdef BL_USE_MPI
    MPI_Finalize();
#endif

    return result;
}