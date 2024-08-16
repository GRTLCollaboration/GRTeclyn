/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef STATEVARIABLES_HPP
#define STATEVARIABLES_HPP

#include "ArrayTools.hpp"
#include "CCZ4Variables.hpp"

/// This enum gives the index of every variable stored in the grid
enum
{
    // Note that it is important that the first enum value is set to 1 more than
    // the last CCZ4 var enum
    NUM_VARS = NUM_CCZ4_VARS,
};

namespace StateVariables
{
static const amrex::Vector<std::string> names = CCZ4Variables::names;
} // namespace StateVariables

#endif /* STATEVARIABLES_HPP */
