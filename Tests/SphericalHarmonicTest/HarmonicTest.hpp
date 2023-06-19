/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef HARMONICTEST_HPP_
#define HARMONICTEST_HPP_

// AMReX includes
#include "AMReX_Array.H"

// #include "BoxLoops.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "Tensor.hpp"
// #include "UserVariables.hpp" //This files needs NUM_VARS - total number of
// components #include "simd.hpp"

class HarmonicTest
{
  public:
    HarmonicTest(std::array<double, AMREX_SPACEDIM> a_center, double a_dx)
        : m_dx(a_dx), m_center(a_center)
    {
    }

    template <class data_t>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
    compute(int i, int j, int k,
            const amrex::CellData<data_t> &current_cell) const;

  private:
    double m_dx;
    std::array<double, AMREX_SPACEDIM> m_center;

    template <class data_t>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE data_t
    compute_harmonic(Coordinates<data_t> coords) const;
};

#include "HarmonicTest.impl.hpp"

#endif /* HARMONICTEST_HPP_ */
