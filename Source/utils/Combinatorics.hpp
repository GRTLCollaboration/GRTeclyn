/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COMBINATORICS_HPP_
#define COMBINATORICS_HPP_

#include <AMReX_BLassert.H>

namespace Combinatorics
{
// Calculate factorials
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double factorial(int n)
{
    double out = 1.0;
    for (int i = 1; i <= n; i++)
    {
        out *= i;
    }
    return out;
}

// Calculate combinatorics
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double n_choose_r(int n, int r)
{
    AMREX_ASSERT((r >= 0) && (n >= r)); // sense check values

    double out = factorial(n) / (factorial(r) * factorial(n - r));
    return out;
}
} // namespace Combinatorics

#endif /* COMBINATORICS_HPP_ */
