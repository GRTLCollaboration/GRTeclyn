/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef MPILAYOUTPARTICLE_HPP_
#define MPILAYOUTPARTICLE_HPP_

#include <vector>

// This class sets up the data structure for MPI communication in MPIContext.
// This is the "plumbing"/"abstract" part of the code: the actual MPI
// communication logic is set-up in MPIContext. If it is not very evident, why
// this is needed, please refer to the explanation inside the MPIContext class.
// This class will keep track of the simplest information relating to each rank,
// i.e. how many points it stores and where each rank's chunk starts in a flat
// buffer.

class MPILayoutParticle
{

  public:
    // Getters
    inline int count(int rank) const;
    inline int totalCount() const;
    inline int displ(int rank) const;

    // Setters
    inline void setCount(int rank, int count);
    inline void incrementCount(int rank);
    inline void clearCounts();

  private:
    friend class MPIContextParticle;

    MPILayoutParticle(int num_process);

    int m_num_process;         // number of processes/ranks
    std::vector<int> m_counts; // how many items/points per rank
    mutable std::vector<int>
        m_displs; // starting offsets in the flat buffer per rank

    mutable int m_total_count{}; // sum of counts over all ranks
    mutable bool m_dirty{false}; // flag to recompute counts/displacements

    inline void updateDirty() const;
    inline int *countsPtr();
    inline int *displsPtr();
};

#include "MPILayoutParticle.impl.hpp"

#endif /* MPILAYOUTPARTICLE_HPP_ */