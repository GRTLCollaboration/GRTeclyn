/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef STATEVARIABLES_HPP
#define STATEVARIABLES_HPP

#include <array>
#include <string>

// assign enum to each variable
enum
{
    c_phi_Re,
    c_phi_Im,

    NUM_VARS
};

namespace StateVariables
{
static const amrex::Vector<std::string> names = {"phi_Re", "phi_Im"};
}

#endif /* STATEVARIABLES_HPP */
