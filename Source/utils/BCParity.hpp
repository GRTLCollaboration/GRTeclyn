/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef BCPARITY_HPP_
#define BCPARITY_HPP_

//! A scoped enumerator for the parity of a variable wrt reflective boundaries
enum class BCParity
{
    even,
    odd_x,
    odd_y,
    odd_z,
    odd_xy,
    odd_yz,
    odd_xz,
    odd_xyz,
    undefined
};

#endif