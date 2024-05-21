/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef USERVARIABLES_HPP
#define USERVARIABLES_HPP

#include "ArrayTools.hpp"
#include "BCParity.hpp"
#include "CCZ4UserVariables.hpp"

#include <array>
#include <string>

enum
{
    c_phi = NUM_CCZ4_VARS,
    c_Pi,

    NUM_VARS
};

namespace UserVariables
{
static const std::array<std::string, NUM_VARS - NUM_CCZ4_VARS>
    user_variable_names = {"phi", "Pi"};

static const std::array<std::string, NUM_VARS> variable_names =
    ArrayTools::concatenate(ccz4_variable_names, user_variable_names);

static const std::array<BCParity, NUM_VARS - NUM_CCZ4_VARS>
    user_variable_parities = {BCParity::even, BCParity::even};

static const std::array<BCParity, NUM_VARS> variable_parities =
    ArrayTools::concatenate(ccz4_variable_parities, user_variable_parities);
} // namespace UserVariables

#endif /* USERVARIABLES_HPP */
