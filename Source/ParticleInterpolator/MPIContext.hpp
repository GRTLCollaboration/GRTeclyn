/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef MPICONTEXT_HPP_
#define MPICONTEXT_HPP_

#include "MPILayout.hpp"
#include <AMReX_ParallelDescriptor.H>
#include <vector>

// This class sets up the MPI communication between the different ranks.
// In particular, this is necessary because our interpolation routine allows for
// a flexible logic, where (i) different ranks can work on (i.e. interpolate)
// different points (ii) different ranks can own different query points This
// means ownership and interpolation of the point may not be on the same rank.
// For example: rank 0 wants to interpolate at point A and rank 1 at point B.
// But: the query of point A is owned by rank 1 and point B by rank 0! (notice
// the swap here) So, when rank 0 "queries" point B, the interpolated value
// needs to be sent back from rank 1; and similarly for rank 1 querying point A.

class MPIContext
{
  public:
    MPIContext();

    // Getters
    inline int queryCount(int rank);
    inline int totalQueryCount();
    inline int answerCount(int rank);
    inline int totalAnswerCount();
    inline int queryDispl(int rank);
    inline int answerDispl(int rank);

    // Setters
    inline void setQueryCount(int rank, int count);
    inline void incrementQueryCount(int rank);
    inline void clearQueryCounts();

    void exchangeLayout();

#ifdef AMREX_USE_MPI
    // MPI asynchronous comms
    inline void asyncBegin();
    inline void asyncExchangeQuery(void *sendbuf, void *recvbuf,
                                   MPI_Datatype type);
    inline void asyncExchangeAnswer(void *sendbuf, void *recvbuf,
                                    MPI_Datatype type);
    inline void asyncEnd();
#endif

    // MPI utils
    static int comm_size();
    static int comm_rank();

  private:
    MPILayout m_query;  // how many things I want to send to each rank
    MPILayout m_answer; // how many things I want to receive from each rank

    bool m_async_active{false};
#ifdef AMREX_USE_MPI
    std::vector<MPI_Request> m_mpi_requests;
#endif
};

inline int MPIContext::comm_size()
{
    return amrex::ParallelDescriptor::NProcs();
}

inline int MPIContext::comm_rank()
{
    return amrex::ParallelDescriptor::MyProc();
}

#include "MPIContext.impl.hpp"

#endif /* MPICONTEXT_HPP_ */
