/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(POLYNOMIALTEST_HPP_)
#error "This file should only be included through PolynomialTest.hpp"
#endif

#ifndef POLYNOMIALTEST_IMPL_HPP_
#define POLYNOMIALTEST_IMPL_HPP_

// This is not pretty and probably redudant, but it works for now.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
PolynomialTest::compute(int i, int j, int k) const
{
    Coordinates coords{
        amrex::IntVect{i, j, k},
        m_dx, m_center
    };

    double poly = compute_polynomial(coords.x, coords.y, coords.z);

    return poly;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
PolynomialTest::compute_polynomial(double x, double y, double z) const
{
    double out = 42. + x * x + y * y * z * z;

    return out;
}

#endif /* POLYNOMIALTEST_IMPL_HPP_ */
