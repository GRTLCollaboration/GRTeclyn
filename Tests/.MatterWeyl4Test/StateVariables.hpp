/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */
#ifndef STATEVARIABLES_HPP
#define STATEVARIABLES_HPP

#include "ArrayTools.hpp"
#include "CCZ4Variables.hpp"

// TODO: This file can be auto-generated from a list of variable names
// Also, we should probably scope this enum too...
//
enum
{
    c_h     = c_h11,
    c_A     = c_A11,
    c_Gamma = c_Gamma1,
    c_shift = c_shift1,
    c_B     = c_B1,

    c_phi = NUM_CCZ4_VARS,
    c_Pi,

    c_Rho,

    c_chi2,

    c_Weyl4_Re,
    c_Weyl4_Im,

    NUM_VARS
};

namespace StateVariables
{
static const amrex::Vector<std::string> additional_names = {

    "phi",      "Pi",

    "rho",

    "chi2",

    "Weyl4_Re", "Weyl4_Im"};

static const amrex::Vector<std::string> names =
    ArrayTools::concatenate(CCZ4Variables::names, additional_names);
} // namespace StateVariables

#endif /* STATEVARIABLES_HPP */
