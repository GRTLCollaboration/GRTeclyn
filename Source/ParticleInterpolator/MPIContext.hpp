/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef MPICONTEXT_HPP_
#define MPICONTEXT_HPP_

#include "MPILayout.hpp"
#include <AMReX_ParallelDescriptor.H>
#include <vector>

// This class sets up the MPI communication between the 'answering' and
// 'quering' ranks. In particular, this is necessary because our interpolation
// routine allows for a flexible logic, where (i) 'answering' ranks can work on
// (i.e. interpolate) different points (ii) 'quering' ranks can ask for some
// different sets of points. This means that asking for a point or interpolating
// a point may happen on DIFFERENT ranks. For example: rank 0 wants to
// interpolate at point A and rank 1 at point B. But: the query of point A
// happens on rank 1 and of point B on rank 0! (notice the swap here) So, when
// rank 0 "queries" point B, the interpolated value needs to be sent back from
// rank 1; and similarly for point A.

// Our philosophy here is therefore -- answering ranks are sending stuff whilst
// quering ranks are receiving stuff.

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
    inline void setAnswerCount(int rank, int count);
    inline void incrementAnswerCount(int rank);
    inline void clearAnswerCounts();

    void exchangeLayout();

#ifdef AMREX_USE_MPI
    // MPI asynchronous comms
    inline void asyncBegin();
    inline void asyncExchangeAnswer(void *sendbuf, void *recvbuf,
                                    MPI_Datatype type);
    inline void asyncEnd();
#endif

    // MPI utils
    static int comm_size();
    static int comm_rank();

  private:
    MPILayout m_query;  // things this ranks wants to receive from each
                        // other rank
    MPILayout m_answer; // things this rank wants to send to each other rank

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