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
    c_A,
    c_B,

    NUM_VARS
};

namespace StateVariables
{
static const std::array<std::string, NUM_VARS> variable_names = {"A", "B"};
}

#endif /* STATEVARIABLES_HPP */
