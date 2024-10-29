/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef STATEVARIABLES_HPP
#define STATEVARIABLES_HPP

#include "ArrayTools.hpp"
#include "BCParity.hpp"

// assign an enum to each variable
enum
{
    // Note that it is important that the first enum value is set to 1 more than

    c_phi,
    c_Pi,

    NUM_VARS = 2,
};

namespace StateVariables
{
static const amrex::Vector<std::string> names{"phi", "Pi"};

static const std::array<BCParity, NUM_VARS> parities = {BCParity::even,
                                                        BCParity::even};
// The example parameter file uses periodic boundary conditions, in which the
// parities aren't used but must be defined. However you could use reflective
// boundary conditions with the above parity definitions.

} // namespace StateVariables
// #include "UserVariables.inc.hpp"

#endif /* STATEVARIABLES_HPP */
