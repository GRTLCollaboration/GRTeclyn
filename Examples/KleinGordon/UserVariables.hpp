/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef USERVARIABLES_HPP
#define USERVARIABLES_HPP

#include "ArrayTools.hpp"
#include "CCZ4UserVariables.hpp"
#include "DiagnosticVariables.hpp"

// assign an enum to each variable                                                                   
enum
{
    // Note that it is important that the first enum value is set to 1 more than                     

    c_phi, 
    c_Pi,   

    NUM_VARS,
};

namespace UserVariables
{
static const std::array<std::string, NUM_VARS>
    variable_names = {"phi", "Pi"};

} // namespace UserVariables     
#include "UserVariables.inc.hpp"

#endif /* USERVARIABLES_HPP */
