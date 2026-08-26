/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SPHERICALHARMONICS_HPP_
#define SPHERICALHARMONICS_HPP_

#include "Combinatorics.hpp"
#include "Coordinates.hpp"
#include "Tensor.hpp"

using namespace amrex::literals;

// Functions for the spin weighted spherical harmonics
// See paper arXiv:gr-qc/0610128 eqn 40
namespace SphericalHarmonics
{
struct Y_lm_t
{
    amrex::Real Real;
    amrex::Real Im;
    amrex::Real magnitude;
};

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
// Calculates the spin weight es, el, em spherical harmonic
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Y_lm_t
spin_Y_lm(const amrex::Real x, const amrex::Real y, const amrex::Real z,
          const int es, const int el, const int em)
{

    AMREX_ASSERT((el >= 0) && (el >= std::abs(em)));

    Y_lm_t Y_lm{};

    // calculate useful position quantities
    amrex::Real r     = sqrt(x * x + y * y + z * z);
    r                 = std::max(r, 1.0e-6_rt);
    amrex::Real theta = acos(z / r);
    amrex::Real phi   = atan2(y, x);

    using namespace Combinatorics;
    amrex::Real coefficient =
        pow(-1.0, es) * sqrt((2.0 * el + 1.0) / (4.0 * M_PI));
    coefficient *= sqrt(factorial(el + em) * factorial(el - em) /
                        factorial(el + es) / factorial(el - es));

    amrex::Real sum = 0.0;
    int lower_limit = em + es > 0 ? em + es : 0;
    int upper_limit = el + em < el + es ? el + em : el + es;

    for (int i = lower_limit; i <= upper_limit; i++)
    {
        amrex::Real temp =
            n_choose_r(el + es, i) * n_choose_r(el - es, i - es - em);
        sum += temp * pow(-1.0, i) *
               pow(cos(theta / 2.0), 2 * (el - i) + es + em) *
               pow(sin(theta / 2.0), 2 * i - em - es);
    }

    Y_lm.Real      = coefficient * sum * cos(em * phi);
    Y_lm.Im        = coefficient * sum * sin(em * phi);
    Y_lm.magnitude = sqrt(Y_lm.Real * Y_lm.Real + Y_lm.Im * Y_lm.Im);

    return Y_lm;
}
// NOLINTEND(bugprone-easily-swappable-parameters)

} // namespace SphericalHarmonics

#endif /* SPHERICALHARMONICS_HPP_ */
