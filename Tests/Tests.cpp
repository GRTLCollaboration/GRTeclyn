#include <iostream>
// Catch2 header
#include "catch_amalgamated.hpp"

// Common test headers
#include "TestsArgs.hpp"

#include "AMReX_REAL.H"
#include "AMReX_ccse-mpi.H"

int main(int argc, char *argv[])
{
    Tests::g_args.set(argc, argv);
#ifdef BL_USE_MPI
    // We can only initialize and finalize MPI once so do it here
    MPI_Init(&argc, &argv);
#endif

    constexpr int cout_precision          = 17;
    Catch::StringMaker<double>::precision = cout_precision;
    Catch::StringMaker<float>::precision  = cout_precision;

    int result = Catch::Session().run(argc, argv);

#ifdef BL_USE_MPI
    MPI_Finalize();
#endif

    return result;
}