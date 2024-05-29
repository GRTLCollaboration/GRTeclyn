/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef USERVARIABLES_HPP
#define USERVARIABLES_HPP

#include "EmptyDiagnosticVariables.hpp"
#include <array>
#include <string>

// assign enum to each variable
enum
{
    c_A,
    c_B,

    NUM_VARS
};

namespace UserVariables
{
static const std::array<std::string, NUM_VARS> variable_names = {"A", "B"};
}

#include "UserVariables.inc.hpp"

#endif /* USERVARIABLES_HPP */
