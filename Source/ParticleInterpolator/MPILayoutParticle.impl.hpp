/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef MPILAYOUTPARTICLE_IMPL_HPP_
#define MPILAYOUTPARTICLE_IMPL_HPP_

inline MPILayoutParticle::MPILayoutParticle(int num_process)
    : m_num_process(num_process), m_counts(m_num_process, 0),
      m_displs(m_num_process, 0)
{
}

// return the number of points for a rank
inline int MPILayoutParticle::count(int rank) const { return m_counts[rank]; }

inline int MPILayoutParticle::totalCount() const
{
    if (m_dirty)
    {
        updateDirty();
    }
    return m_total_count;
}

// return the offset in the flat buffer for a rank
inline int MPILayoutParticle::displ(int rank) const
{
    if (m_dirty)
    {
        updateDirty();
    }
    return m_displs[rank];
}

inline void MPILayoutParticle::setCount(int rank, int count)
{
    AMREX_ASSERT(rank < m_num_process && count >= 0);
    m_counts[rank] = count;
    m_dirty        = true;
}

// increment by one the number of points for a rank
inline void MPILayoutParticle::incrementCount(int rank)
{
    AMREX_ASSERT(rank < m_num_process);
    ++m_counts[rank];
    m_dirty = true;
}

inline void MPILayoutParticle::clearCounts()
{
    m_counts.assign(m_num_process, 0);
    m_dirty = true;
}

// This function essential divvides the flat buffer into per rank segments,
// which have their starting points and displacements
inline void MPILayoutParticle::updateDirty() const
{
    m_total_count = m_counts[0];
    for (int i = 1; i < m_num_process; ++i)
    {
        m_total_count += m_counts[i];
        m_displs[i]    = m_displs[i - 1] + m_counts[i - 1];
    }
    m_dirty = false;
}

inline int *MPILayoutParticle::countsPtr() { return m_counts.data(); }

inline int *MPILayoutParticle::displsPtr()
{
    if (m_dirty)
    {
        updateDirty();
    }
    return m_displs.data();
}

#endif /* MPILAYOUTPARTICLE_IMPL_HPP_ */
