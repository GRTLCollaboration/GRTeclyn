/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef CCZ4GEOMETRYUNITTEST_HPP_
#define CCZ4GEOMETRYUNITTEST_HPP_

#include "DimensionDefinitions.hpp"

// We test chris.ULL (rank 3) + h_UU and RicciZ_.LL (2 x rank 2)
// + chris.contracted (rank 1) + RicciZ.scalar (scalar)
constexpr int NUM_GEOMETRY_TEST_VARS = GR_SPACEDIM * GR_SPACEDIM * GR_SPACEDIM +
                                       2 * GR_SPACEDIM * GR_SPACEDIM +
                                       GR_SPACEDIM + 1;

void run_ccz4_geometry_unit_tests();

#endif /* CCZ4GEOMETRYUNITTEST_HPP_ */