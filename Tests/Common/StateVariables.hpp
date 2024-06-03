/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef STATEVARIABLES_HPP
#define STATEVARIABLES_HPP

#include "ArrayTools.hpp"
#include "BCParity.hpp"
#include "CCZ4Variables.hpp"

#include <array>
#include <string>

enum
{
    c_phi = NUM_CCZ4_VARS,
    c_Pi,

    NUM_VARS
};

namespace StateVariables
{
static const amrex::Vector<std::string> additional_names = {"phi", "Pi"};

static const amrex::Vector<std::string> names =
    ArrayTools::concatenate(CCZ4Variables::names, additional_names);

static const std::array<BCParity, NUM_VARS - NUM_CCZ4_VARS>
    user_variable_parities = {BCParity::even, BCParity::even};

static const std::array<BCParity, NUM_VARS> parities =
    ArrayTools::concatenate(CCZ4Variables::parities, user_variable_parities);
} // namespace StateVariables

#endif /* STATEVARIABLES_HPP */
