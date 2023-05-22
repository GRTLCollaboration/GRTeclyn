/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef USERVARIABLES_HPP
#define USERVARIABLES_HPP

#include "ArrayTools.hpp"
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
} // namespace UserVariables

#endif /* USERVARIABLES_HPP */
