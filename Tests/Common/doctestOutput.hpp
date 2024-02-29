#ifndef DOCTESTOUTPUT_HPP_
#define DOCTESTOUTPUT_HPP_

// Doctest header
#include "doctest.h"

// AMReX includes
#include "AMReX_ccse-mpi.H"

// System includes
#include <fstream>
#include <iostream>

#ifdef BL_USE_MPI
namespace doctest
{
// Hide output from non-zero ranks when building with MPI
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
} // namespace doctest
#endif

#endif /* DOCTESTOUTPUT_HPP_ */