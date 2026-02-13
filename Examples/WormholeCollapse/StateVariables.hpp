#ifndef STATEVARIABLES_HPP
#define STATEVARIABLES_HPP

#include "CCZ4StateVariables.hpp"

// We use the standard CCZ4 variables (Chi, Metric, Gamma, K, A, Theta, Lapse,
// Shift, B)
enum
{
    NUM_VARS = NUM_CCZ4_VARS,
};

namespace StateVariables
{
static const amrex::Vector<std::string> names = CCZ4StateVariables::names;
static const std::array<BCParity, NUM_VARS> parities =
    CCZ4StateVariables::parities;
} // namespace StateVariables

#endif /* STATEVARIABLES_HPP */