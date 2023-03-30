/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INTERPOLATIONALGORITHM_HPP_
#define INTERPOLATIONALGORITHM_HPP_

#include <AMReX_FArrayBox.H>

// Other includes
#include <array>

class InterpolationAlgorithm
{
};

class NearestNeighbour : public InterpolationAlgorithm
{
  public:
    static inline double
    interpPoint(const std::array<double, AMREX_SPACEDIM> & /*gridCoord*/,
                const amrex::FArrayBox & /*fab*/, int /*comps*/,
                const amrex::IntVect & /*nearest*/)
    {
        amrex::Abort("xxxxx interpPoint todo");
        return 0.; // xxxxx fab.get(nearest, comps);
    }
};

#endif /* INTERPOLATIONALGORITHM_HPP_ */
