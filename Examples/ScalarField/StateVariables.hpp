/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef STATEVARIABLES_HPP
#define STATEVARIABLES_HPP

#include "ArrayTools.hpp"
#include "BCParity.hpp"

#include "CCZ4Variables.hpp"

// assign an enum to each variable
enum
{
    // Note that it is important that the first enum value is set to 1 more than
    // the last CCZ4 var enum
    c_phi = NUM_CCZ4_VARS, // matter field added
    c_Pi,                  //(minus) conjugate momentum

    NUM_VARS
};

namespace StateVariables
{
static const std::array<std::string, NUM_VARS - NUM_CCZ4_VARS>
    additional_names = {"phi", "Pi"};

static const std::array<std::string, NUM_VARS> names =
    ArrayTools::concatenate(CCZ4Variables::names, additional_names);

static const std::array<BCParity, NUM_VARS - NUM_CCZ4_VARS>
    additional_parities = {BCParity::even, BCParity::even};

static const std::array<BCParity, NUM_VARS> parities =
    ArrayTools::concatenate(CCZ4Variables::parities, additional_parities);

} // namespace StateVariables

#endif /* STATEVARIABLES_HPP */
