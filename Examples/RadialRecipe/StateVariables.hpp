#ifndef STATEVARIABLES_HPP
#define STATEVARIABLES_HPP

#include "ArrayTools.hpp"
#include "CCZ4StateVariables.hpp"
#include "GRTresnaScalarLayout.hpp"

enum
{
    c_phi  = NUM_CCZ4_VARS,
    c_Pi,
    // Complex scalar reuses the first lump slots (phi_lump0 / Pi_lump0).
    c_phi2 = NUM_CCZ4_VARS + 2,
    c_Pi2  = NUM_CCZ4_VARS + 3,

    NUM_VARS = NUM_CCZ4_VARS + 2 +
               2 * GRTRESNA_MAX_INDEPENDENT_SCALARS
};

namespace StateVariables
{
static const amrex::Vector<std::string> additional_names = {
    "phi", "Pi",       "phi_lump0", "Pi_lump0", "phi_lump1", "Pi_lump1",
    "phi_lump2", "Pi_lump2", "phi_lump3", "Pi_lump3", "phi_lump4", "Pi_lump4"};
static const amrex::Vector<std::string> names =
    ArrayTools::concatenate(CCZ4StateVariables::names, additional_names);

static const std::array<BCParity, 12> additional_parities = {
    BCParity::even, BCParity::even, BCParity::even, BCParity::even,
    BCParity::even, BCParity::even, BCParity::even, BCParity::even,
    BCParity::even, BCParity::even, BCParity::even, BCParity::even};
static const std::array<BCParity, NUM_VARS> parities =
    ArrayTools::concatenate(CCZ4StateVariables::parities, additional_parities);
} // namespace StateVariables

#endif /* STATEVARIABLES_HPP */
