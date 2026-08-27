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

static constexpr int NUM_MATTER_VARS = NUM_VARS - NUM_CCZ4_VARS;

static const std::array<BCParity, NUM_MATTER_VARS> scalar_field_parities = {
    BCParity::even, BCParity::even};

static const std::array<BCParity, NUM_VARS> parities = ArrayTools::concatenate(
    CCZ4StateVariables::parities, scalar_field_parities);

static const std::array<amrex::Real, NUM_MATTER_VARS>
    scalar_field_asymptotic_values = {0.0, 0.0};

static const std::array<amrex::Real, NUM_VARS> asymptotic_values =
    ArrayTools::concatenate(CCZ4StateVariables::asymptotic_values,
                            scalar_field_asymptotic_values);
} // namespace StateVariables

#endif /* STATEVARIABLES_HPP */
