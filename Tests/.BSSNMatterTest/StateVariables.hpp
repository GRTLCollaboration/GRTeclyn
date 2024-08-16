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
    c_h     = c_h11,
    c_A     = c_A11,
    c_Gamma = c_Gamma1,
    c_shift = c_shift1,
    c_B     = c_B1,
    // Note that it is important that the first enum value is set to 1 more than
    // the last CCZ4 var enum
    c_phi = NUM_CCZ4_VARS,
    c_Pi,

    c_chi2,

    c_Ham,

    c_Mom,
    c_Mom1 = c_Mom,
    c_Mom2,
    c_Mom3,

    NUM_VARS
};

namespace StateVariables
{
static const amrex::Vector<std::string> additional_names = {"phi",  "Pi",

                                                            "chi2",

                                                            "Ham",

                                                            "Mom1", "Mom2",
                                                            "Mom3"};

static const amrex::Vector<std::string> names =
    ArrayTools::concatenate(CCZ4Variables::names, additional_names);
} // namespace StateVariables

#endif /* STATEVARIABLES_HPP */
