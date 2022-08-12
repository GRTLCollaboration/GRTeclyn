/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef QUINTICCONVOLUTION_HPP_
#define QUINTICCONVOLUTION_HPP_

#include "InterpSource.hpp"
#include <utility>

class QuinticConvolution
{
    const InterpSource &m_source;
    bool m_verbosity;

    std::vector<IntVect> m_interp_points;
    std::vector<double> m_interp_weights;

  public:
    QuinticConvolution(const InterpSource &source, bool verbosity = false);

    void setup(const std::array<int, AMREX_SPACEDIM> &deriv,
               const std::array<double, AMREX_SPACEDIM> &dx,
               const std::array<double, AMREX_SPACEDIM> &evalCoord,
               const IntVect &nearest);
    double interpData(const FArrayBox &fab, int comp);

    const static string TAG;
};

#include "QuinticConvolution.impl.hpp"

#endif /* QUINTICCONVOLUTION_HPP_ */
