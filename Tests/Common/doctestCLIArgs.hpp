/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef DOCTESTCLIARGS_HPP_
#define DOCTESTCLIARGS_HPP_

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

namespace doctest
{

//! A class to remove prefixed doctest options from the command line args
class CLIArgs
{
    std::vector<char *> m_args_vec;

  public:
    void set(char **a_argv)
    {
        for (; *a_argv != nullptr; ++a_argv)
        {
            if (std::string(*a_argv).find("-dt-") == std::string::npos)
            {
                m_args_vec.push_back(*a_argv);
            }
        }
        m_args_vec.push_back(nullptr);
    }

    int argc() { return static_cast<int>(m_args_vec.size()) - 1; }
    char **argv() { return m_args_vec.data(); } // Note: non-const char **:

    //! Insert an argument at a specific index (e.g. to add a parameter file)
    void insert(int a_arg_index, char *a_arg)
    {
        // negative indexes start at the end but leave the terminating nullptr
        auto pos_it =
            (a_arg_index > 0) ? m_args_vec.begin() : m_args_vec.end() - 1;

        // make sure we're not going past the end of the vector
        a_arg_index = std::clamp(a_arg_index, -argc(), argc());

        pos_it += a_arg_index;

        m_args_vec.insert(pos_it, a_arg);
    }
};

// Unfortunately the following has to be global and non-const
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
extern CLIArgs cli_args;
}; // namespace doctest

#endif /* DOCTESTCLIARGS_HPP_ */