/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef POLYNOMIALTEST_HPP_
#define POLYNOMIALTEST_HPP_

// AMReX includes
#include "AMReX_Array.H"

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "Tensor.hpp"

class PolynomialTest
{
  public:
    PolynomialTest(std::array<double, AMREX_SPACEDIM> a_center, double a_dx)
        : m_dx(a_dx), m_center(a_center)
    {
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double compute(int i, int j,
                                                            int k) const;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
    compute_polynomial(double x, double y, double z) const;

  private:
    double m_dx;
    std::array<double, AMREX_SPACEDIM> m_center;
};

#include "PolynomialTest.impl.hpp"

#endif /* POLYNOMIALTEST_HPP_ */
