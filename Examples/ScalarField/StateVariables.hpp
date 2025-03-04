/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef STATEVARIABLES_HPP
#define STATEVARIABLES_HPP

#include "ArrayTools.hpp"
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
static const amrex::Vector<std::string> scalar_field_names = {"phi", "Pi"};

static const amrex::Vector<std::string> names =
    ArrayTools::concatenate(CCZ4Variables::names, scalar_field_names);

static const std::array<BCParity, NUM_VARS - NUM_CCZ4_VARS>
    scalar_field_parities = {BCParity::even, BCParity::even};

static const std::array<BCParity, NUM_VARS> parities =
    ArrayTools::concatenate(CCZ4Variables::parities, scalar_field_parities);
} // namespace StateVariables

#endif /* STATEVARIABLES_HPP */
