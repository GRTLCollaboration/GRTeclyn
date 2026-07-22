/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef CCZ4GEOMETRYUNITTEST_HPP_
#define CCZ4GEOMETRYUNITTEST_HPP_

#include "DimensionDefinitions.hpp"

// We test chris.ULL (rank 3) + chris_phys.ULL (rank3) + h_UU and RicciZ_.LL (2
// x rank 2)
// + chris.contracted (rank 1) + RicciZ.scalar (scalar)
// + RicciZ_general.LL (2 x rank 2) + RicciZ_general.scalar (scalar)
constexpr int NUM_GEOMETRY_TEST_VARS =
    2 * AMREX_SPACEDIM * AMREX_SPACEDIM * AMREX_SPACEDIM +
    2 * AMREX_SPACEDIM * AMREX_SPACEDIM + AMREX_SPACEDIM + 1 +
    AMREX_SPACEDIM * AMREX_SPACEDIM + 1;

void run_ccz4_geometry_unit_tests();

#endif /* CCZ4GEOMETRYUNITTEST_HPP_ */
