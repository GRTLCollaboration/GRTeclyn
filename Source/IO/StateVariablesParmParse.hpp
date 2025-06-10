/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef STATEVARIABLESPARMPARSE_HPP
#define STATEVARIABLESPARMPARSE_HPP

#include "StateVariables.hpp"

// Other includes
#include "GRParmParse.hpp"
#include "VariableType.hpp"
#include <algorithm>
#include <array>
#include <string>

//! Some functions to help loading state variables from a ParmParse
namespace StateVariablesParmParse
{
/// Takes a string and returns the variable enum number if the string
/// matches one of those in StateVariables::names, or returns -1
/// otherwise
static int name_to_enum(const std::string &a_var_name)
{
    const auto var_name_it =
        std::find(StateVariables::names.cbegin(), StateVariables::names.cend(),
                  a_var_name);

    auto var = std::distance(StateVariables::names.begin(), var_name_it);
    if (var != NUM_VARS)
    {
        return static_cast<int>(var);
    }
    return -1;
}

// where one has read in a subset of variables with some feature
// this reads in a set of associated values and assigns it into a full
// array of all NUM_VARS vars (setting other values to a default value)
template <class T>
void load_values_to_array(GRParmParse &pp, const char *a_values_vector_string,
                          const std::vector<int> &a_vars_vector,
                          std::array<T, NUM_VARS> &a_values_array,
                          const T a_default_value)
{
    // how many values do I need to get?
    auto num_values = a_vars_vector.size();
    // make a container for them, and load
    std::vector<T> vars_values(num_values, a_default_value);
    pp.load(a_values_vector_string, vars_values, num_values, vars_values);

    // populate the values_array for the NUM_VARS values with those read in
    a_values_array.fill(a_default_value);
    for (int i = 0; i < num_values; i++)
    {
        int icomp             = a_vars_vector[i];
        a_values_array[icomp] = vars_values[i];
    }
}

// function to create a vector of enums of vars by reading in their
// names as strings from the params file and converting it to the enums
inline void load_vars_to_vector(GRParmParse &pp,
                                const char *a_vars_vector_string,
                                std::vector<int> &a_vars_vector)
{
    int num_values = pp.countval(a_vars_vector_string);
    if (num_values > 0)
    {
        std::vector<std::string> var_names(num_values, "");
        pp.load(a_vars_vector_string, var_names, num_values, var_names);
        for (const std::string &var_name : var_names)
        {
            int var = name_to_enum(var_name);
            if (var < 0)
            {
                amrex::Print()
                    << "Variable with name " << var_name << " not found.\n";
            }
            else
            {
                a_vars_vector.push_back(var);
            }
        }
    }
}

} // namespace StateVariablesParmParse

#endif /* STATEVARIABLESPARMPARSE_HPP */
