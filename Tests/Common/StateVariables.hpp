/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef STATEVARIABLES_HPP
#define STATEVARIABLES_HPP

#include "ArrayTools.hpp"
#include "BCParity.hpp"
#include "CCZ4StateVariables.hpp"

#include <array>
#include <string>

enum
{
    c_phi = NUM_CCZ4_VARS,
    c_Pi,
    c_phi_Re, // Re(Y_{20})
    c_phi_Im, // Im(Y_{20}) = 0 in here
    c_polystate,

    NUM_VARS
};

namespace StateVariables
{
static const amrex::Vector<std::string> additional_names = {
    "phi", "Pi", "phi_Re", "phi_Im", "polystate"};

static const amrex::Vector<std::string> names =
    ArrayTools::concatenate(CCZ4StateVariables::names, additional_names);

static const std::array<BCParity, NUM_VARS - NUM_CCZ4_VARS>
    user_variable_parities = {BCParity::even, BCParity::even, BCParity::even,
                              BCParity::even, BCParity::odd_z};

static const std::array<BCParity, NUM_VARS> parities = ArrayTools::concatenate(
    CCZ4StateVariables::parities, user_variable_parities);
} // namespace StateVariables

#endif /* STATEVARIABLES_HPP */
