/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef USERVARIABLES_HPP
#define USERVARIABLES_HPP

#include <array>
#include <string>

// assign enum to each variable
enum
{
    c_phi_Re,
    c_phi_Im,

    NUM_VARS
};

namespace UserVariables
{
static const std::array<std::string, NUM_VARS> variable_names = {"phi_Re",
                                                                 "phi_Im"};
}

#include "UserVariables.inc.hpp"

#endif /* USERVARIABLES_HPP */
