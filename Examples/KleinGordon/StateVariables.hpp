/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef STATEVARIABLES_HPP
#define STATEVARIABLES_HPP

#include "ArrayTools.hpp"
#include "BCParity.hpp"

// assign an enum to each variable
enum
{
    c_phi,
    c_Pi,

    NUM_VARS,
};

namespace StateVariables
{
static const amrex::Vector<std::string> names{"phi", "Pi"};

static const std::array<BCParity, NUM_VARS> parities = {BCParity::even,
                                                        BCParity::even};
static const std::array<amrex::Real, NUM_VARS> asymptotic_values{};
// The example parameter file uses periodic boundary conditions, in which the
// parities aren't used but must be defined. However you could use reflective
// boundary conditions with the above parity definitions.

} // namespace StateVariables

#endif /* STATEVARIABLES_HPP */
