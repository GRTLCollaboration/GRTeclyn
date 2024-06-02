/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */
#ifndef TESTCASES_HPP_
#define TESTCASES_HPP_

// doctest header
#include "doctest.h"

// AMReX includes
#include <AMReX.H>

// Test cases
#include "CCZ4GeometryUnitTest.hpp"
#include "CCZ4RHSTest.hpp"
#include "ConstraintsTest.hpp"
#include "CoordinateTransformationsTest.hpp"
#include "DerivativeUnitTests.hpp"
#include "MatterCCZ4RHSTest.hpp"
#include "MatterWeyl4Test.hpp"
#include "PositiveChiAndAlphaUnitTest.hpp"
#include "SphericalHarmonicTest.hpp"
#include "Weyl4Test.hpp"

TEST_CASE("CCZ4 Geometry") { run_ccz4_geometry_unit_tests(); }

TEST_CASE("CCZ4 RHS") { run_ccz4_rhs_test(); }

TEST_CASE("Constraints"
#ifndef AMREX_USE_HDF5
          * doctest::skip()
#endif
)
{
    run_constraints_test();
}

TEST_CASE("Coordinate Transformations")
{
    run_coordinate_transformations_test();
}

TEST_CASE("Derivative Unit Tests") { run_derivative_unit_tests(); }

TEST_CASE("Positive Chi and Alpha") { run_positive_chi_and_alpha_unit_test(); }

TEST_CASE("Spherical Harmonics") { run_spherical_harmonic_test(); }

TEST_CASE("Matter Weyl4") { run_matter_weyl4_test(); }

TEST_CASE("Matter CCZ4") { run_matter_ccz4_rhs_test(); }

TEST_CASE("Weyl4"
#ifndef AMREX_USE_HDF5
          * doctest::skip()
#endif
)
{
    run_weyl4_test();
}

#endif /* TESTCASES_HPP_ */
