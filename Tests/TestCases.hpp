/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */
#ifndef TESTCASES_HPP_
#define TESTCASES_HPP_

// doctest header
#include "doctest.h"

// Test cases
#include "CCZ4GeometryUnitTest.hpp"
#include "CCZ4RHSTest.hpp"
#include "CoordinateTransformationsTest.hpp"
#include "DerivativeUnitTests.hpp"
#include "PositiveChiAndAlphaUnitTest.hpp"
#include "SphericalHarmonicTest.hpp"

TEST_CASE("CCZ4 Geometry") { run_ccz4_geometry_unit_tests(); }

TEST_CASE("CCZ4 RHS") { run_ccz4_rhs_test(); }

TEST_CASE("Coordinate Transformations")
{
    run_coordinate_transformations_test();
}

TEST_CASE("Derivative Unit Tests") { run_derivative_unit_tests(); }

TEST_CASE("Positive Chi and Alpha") { run_positive_chi_and_alpha_unit_test(); }

TEST_CASE("Spherical Harmonics") { run_spherical_harmonic_test(); }

#endif /* TESTCASES_HPP_ */
