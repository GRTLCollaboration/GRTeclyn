/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef MPICONTEXTPARTICLE_IMPL_HPP_
#define MPICONTEXTPARTICLE_IMPL_HPP_

inline MPIContextParticle::MPIContextParticle()
    : m_query(comm_size()), m_answer(comm_size())
{
}

inline int MPIContextParticle::queryCount(int rank)
{
    return m_query.count(rank);
}

inline int MPIContextParticle::totalQueryCount()
{
    return m_query.totalCount();
}
inline int MPIContextParticle::answerCount(int rank)
{
    return m_answer.count(rank);
}

inline int MPIContextParticle::totalAnswerCount()
{
    return m_answer.totalCount();
}

inline int MPIContextParticle::queryDispl(int rank)
{
    return m_query.displ(rank);
} // where in the receive buffer the rank's data starts

inline int MPIContextParticle::answerDispl(int rank)
{
    return m_answer.displ(rank);
} // where in the sending buffer the ranks' data starts

inline void MPIContextParticle::setAnswerCount(int rank, int count)
{
    AMREX_ASSERT(!m_async_active);
    m_answer.setCount(rank, count);
}

// how many things I want to send the rank in the arg
inline void MPIContextParticle::incrementAnswerCount(int rank)
{
    AMREX_ASSERT(!m_async_active);
    m_answer.incrementCount(rank);
}

inline void MPIContextParticle::clearAnswerCounts()
{
    AMREX_ASSERT(!m_async_active);
    m_answer.clearCounts();
}

// set up the exchange of points between all the ranks
// m_query is things I want to receive
// m_answer is things I will send
// see also here:
// (https://www.mpich.org/static/docs/v3.1/www3/MPI_Alltoall.html)
inline void MPIContextParticle::exchangeLayout()
{
    AMREX_ASSERT(!m_async_active);
#ifdef AMREX_USE_MPI
    MPI_Alltoall(m_answer.countsPtr(), 1, MPI_INT, m_query.countsPtr(), 1,
                 MPI_INT, amrex::ParallelDescriptor::Communicator());
#else
    *m_query.countsPtr() = *m_answer.countsPtr();
#endif
    m_query.updateDirty();
}

#ifdef AMREX_USE_MPI
inline void MPIContextParticle::asyncBegin()
{
    AMREX_ASSERT(!m_async_active);
    m_async_active = true;
}

inline void MPIContextParticle::asyncExchangeAnswer(void *sendbuf,
                                                    void *recvbuf,
                                                    MPI_Datatype type)
{
    AMREX_ASSERT(m_async_active);
    MPI_Request req = 0;
    m_mpi_requests.push_back(req);

#if MPI_VERSION >= 3
    MPI_Ialltoallv(sendbuf, m_answer.countsPtr(), m_answer.displsPtr(), type,
                   recvbuf, m_query.countsPtr(), m_query.displsPtr(), type,
                   amrex::ParallelDescriptor::Communicator(),
                   &m_mpi_requests.back());
#else
    MPI_Alltoallv(sendbuf, m_answer.countsPtr(), m_answer.displsPtr(), type,
                  recvbuf, m_query.countsPtr(), m_query.displsPtr(), type,
                  amrex::ParallelDescriptor::Communicator());
#endif
}

inline void MPIContextParticle::asyncEnd()
{
    AMREX_ASSERT(m_async_active);
    m_async_active = false;

#if MPI_VERSION >= 3
    MPI_Waitall(
        static_cast<int>(m_mpi_requests.size()), m_mpi_requests.data(),
        MPI_STATUSES_IGNORE); // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
#endif

    m_mpi_requests.clear();
}
#endif /* ifdef AMREX_USE_MPI */

#endif /* MPICONTEXTPARTICLE_IMPL_HPP_ */
