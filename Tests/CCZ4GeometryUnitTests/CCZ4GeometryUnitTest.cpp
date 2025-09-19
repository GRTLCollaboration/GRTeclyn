/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// Doctest header
#include "doctest.h"

// Test header
#include "CCZ4GeometryUnitTest.hpp"

// System includes
// #include <iostream>

// Our includes
#include "CCZ4Geometry.hpp"
#include "DimensionDefinitions.hpp"

template <class data_t> struct vars_t
{
    data_t chi;
    amrex::Array2D<data_t, 0, 3, 0, 3> h;
    amrex::Array1D<data_t, 0, 3> Gamma;
};

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
void run_ccz4_geometry_unit_tests()
{
    vars_t<amrex::Real> vars{};
    vars_t<amrex::Array1D<amrex::Real, 0, 3>> d1{};
    vars_t<amrex::Array2D<amrex::Real, 0, 3, 0, 3>> d2{};
    amrex::Array1D<amrex::Real, 0, 3> Z_over_chi{};

#include "CCZ4GeometryMathematicaValues.hpp" //Including the auto generated file with values

    auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);

    auto chris = TensorAlgebra::compute_christoffel(d1.h, h_UU);

    auto ricciZ =
        CCZ4Geometry::compute_ricci_Z(vars, d1, d2, h_UU, chris, Z_over_chi);

    double test_threshold = 1e-14;

    // Compare
    FOR (i, j)
    {
        INFO("h_UU[" << i << "][" << j << "]");
        CHECK(h_UU(i, j) ==
              doctest::Approx(h_UU_known[i][j]).epsilon(test_threshold));
    }

    FOR (i, j, k)
    {
        INFO("chris.ULL[" << i << "][" << j << "][" << k << "]");
        CHECK(chris.ULL(i, j, k) ==
              doctest::Approx(chris_known[i][j][k]).epsilon(test_threshold));
    }

    FOR (i)
    {
        INFO("chris.contracted[" << i << "]");
        CHECK(
            chris.contracted(i) ==
            doctest::Approx(chris_contracted_known[i]).epsilon(test_threshold));
    }

    FOR (i, j)
    {
        INFO("ricciZ.LL[" << i << "][" << j << "]");
        CHECK(ricciZ.LL(i, j) ==
              doctest::Approx(ricciZ_known[i][j]).epsilon(test_threshold));
    }

    CHECK(ricciZ.scalar ==
          doctest::Approx(ricciZ_scalar_known).epsilon(test_threshold));
}
