/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef CCZ4VARIABLES_HPP
#define CCZ4VARIABLES_HPP

#include <algorithm>
#include <array>
#include <string>

#include "BCParity.hpp"

/// This enum gives the index of the CCZ4 variables on the grid
enum
{
    c_chi,

    c_h11,
    c_h12,
    c_h13,
    c_h22,
    c_h23,
    c_h33,

    c_K,

    c_A11,
    c_A12,
    c_A13,
    c_A22,
    c_A23,
    c_A33,

    c_Theta,

    c_Gamma1,
    c_Gamma2,
    c_Gamma3,

    c_lapse,

    c_shift1,
    c_shift2,
    c_shift3,

    c_B1,
    c_B2,
    c_B3,

    NUM_CCZ4_VARS
};

namespace CCZ4Variables
{
static const amrex::Vector<std::string> names = {
    "chi",

    "h11",    "h12",    "h13",    "h22", "h23", "h33",

    "K",

    "A11",    "A12",    "A13",    "A22", "A23", "A33",

    "Theta",

    "Gamma1", "Gamma2", "Gamma3",

    "lapse",

    "shift1", "shift2", "shift3",

    "B1",     "B2",     "B3",
};

static const std::array<BCParity, NUM_CCZ4_VARS> parities = {
    BCParity::even, // chi

    BCParity::even,   // h11
    BCParity::odd_xy, // h12
    BCParity::odd_xz, // h13
    BCParity::even,   // h22
    BCParity::odd_yz, // h23
    BCParity::even,   // h33

    BCParity::even, // K

    BCParity::even,   // A11
    BCParity::odd_xy, // A12
    BCParity::odd_xz, // A13
    BCParity::even,   // A22
    BCParity::odd_yz, // A23
    BCParity::even,   // A33

    BCParity::even, // Theta

    BCParity::odd_x, // Gamma1
    BCParity::odd_y, // Gamma2
    BCParity::odd_z, // Gamma3

    BCParity::even, // lapse

    BCParity::odd_x, // shift1
    BCParity::odd_y, // shift2
    BCParity::odd_z, // shift3

    BCParity::odd_x, // B1
    BCParity::odd_y, // B2
    BCParity::odd_z, // B3
};
} // namespace CCZ4Variables

#endif /* CCZ4VARIABLES_HPP */
