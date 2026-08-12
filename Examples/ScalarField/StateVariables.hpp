/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef STATEVARIABLES_HPP
#define STATEVARIABLES_HPP

#include "ArrayTools.hpp"
#include "BCParity.hpp"
#include "CCZ4StateVariables.hpp"

/// This enum gives the index of every variable stored in the grid.
enum
{
    c_phi = NUM_CCZ4_VARS,
    c_Pi,

    NUM_VARS
};

namespace StateVariables
{
static const amrex::Vector<std::string> scalar_field_names = {"phi", "Pi"};

static const amrex::Vector<std::string> names =
    ArrayTools::concatenate(CCZ4StateVariables::names, scalar_field_names);

static const std::array<BCParity, 2> scalar_field_parities = {BCParity::even,
                                                              BCParity::even};

static const std::array<BCParity, NUM_VARS> parities = ArrayTools::concatenate(
    CCZ4StateVariables::parities, scalar_field_parities);
} // namespace StateVariables

#endif /* STATEVARIABLES_HPP */
