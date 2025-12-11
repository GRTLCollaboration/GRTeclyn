/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef INTERPOLATIONLAYOUTPARTICLE_HPP_
#define INTERPOLATIONLAYOUTPARTICLE_HPP_

#include <vector>

// This class set-ups a simple storage layour for interpolation queries for
// particles.

template <int num_components> class ParticleInterpolator; // forward declaration

class InterpolationLayoutParticle
{
  private:
    template <int num_components> friend class ParticleInterpolator;

    std::vector<int> rank; // owner rank for each query point
    std::vector<int>
        q_local; // local index of the query point on the owner rank

    InterpolationLayoutParticle(size_t num_points)
        : rank(num_points, -1), q_local(num_points, -1)
    {
    }
};

#endif /* INTERPOLATIONLAYOUTPARTICLE_HPP_ */
