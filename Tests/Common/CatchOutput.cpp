// AMReX includes
#include "AMReX_ccse-mpi.H"

// System includes
#include <fstream>
#include <iostream>

#ifdef BL_USE_MPI
namespace Catch
{
// Catch2 doesn't really know anything about MPI so override the Catch::cout()
// function here so that only rank 0 prints to stdout and the output is not
// garbled.
std::ostream &hide_output_from_non_zero_ranks(std::ostream &a_rank_zero_ostream)
{
    int mpi_initialized = false;
    MPI_Initialized(&mpi_initialized);
    if (!mpi_initialized)
    {
        return a_rank_zero_ostream;
    }
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
    {
        return a_rank_zero_ostream;
    }
    static std::ofstream null_stream;
    if (!null_stream.is_open())
    {
        null_stream.open("/dev/null", std::ofstream::out | std::ofstream::app);
    }
    return null_stream;
}

std::ostream &cout() { return hide_output_from_non_zero_ranks(std::cout); }

std::ostream &cerr() { return std::cerr; }

std::ostream &clog() { return std::clog; }

} // namespace Catch
#endif