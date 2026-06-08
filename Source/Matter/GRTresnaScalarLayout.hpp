#ifndef GRTRESNA_SCALAR_LAYOUT_HPP_
#define GRTRESNA_SCALAR_LAYOUT_HPP_

#include "CCZ4StateVariables.hpp"

//! Maximum independent scalar lumps carried in RadialRecipe state data.
static constexpr int GRTRESNA_MAX_INDEPENDENT_SCALARS = 5;

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr int
c_phi_lump_index(int k)
{
    return NUM_CCZ4_VARS + 2 + 2 * k;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr int
c_Pi_lump_index(int k)
{
    return c_phi_lump_index(k) + 1;
}

#endif /* GRTRESNA_SCALAR_LAYOUT_HPP_ */
