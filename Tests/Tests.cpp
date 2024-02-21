#include <iostream>
// Catch2 header
#include "catch_amalgamated.hpp"

#include "AMReX.H"
#include "AMReX_REAL.H"
#include "AMReX_ccse-mpi.H"

int main(int argc, char *argv[])
{
#ifdef BL_USE_MPI
    // We can only initialize and finalize MPI once so do it here
    MPI_Init(&argc, &argv);
#endif

    Catch::Session catch_session;

    // Default AMReX verbosity to 0 to avoid the "Initialized", "Finalized" and
    // memory usage messages
    amrex::system::verbose = 0;

    auto cli = catch_session.cli() |
               Catch::Clara::Opt(amrex::system::verbose,
                                 "amrex-verbosity")["--amrex-verbosity"](
                   "AMReX Verbosity (default: 0)");

    catch_session.cli(cli);

    constexpr int cout_precision          = 17;
    Catch::StringMaker<double>::precision = cout_precision;
    Catch::StringMaker<float>::precision  = cout_precision;

    int result = catch_session.run(argc, argv);

#ifdef BL_USE_MPI
    MPI_Finalize();
#endif

    return result;
}