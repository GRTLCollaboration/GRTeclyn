#ifndef STATEVARIABLES_HPP
#define STATEVARIABLES_HPP

#include "ArrayTools.hpp"
#include "CCZ4StateVariables.hpp"

enum
{
    c_phi = NUM_CCZ4_VARS,
    c_Pi,

    NUM_VARS
};

namespace StateVariables
{
static const amrex::Vector<std::string> additional_names = {"phi", "Pi"};
static const amrex::Vector<std::string> names =
    ArrayTools::concatenate(CCZ4StateVariables::names, additional_names);

static const std::array<BCParity, 2> additional_parities = {BCParity::even,
                                                            BCParity::even};
static const std::array<BCParity, NUM_VARS> parities =
    ArrayTools::concatenate(CCZ4StateVariables::parities, additional_parities);
} // namespace StateVariables

#endif /* STATEVARIABLES_HPP */
