/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// Catch2 header
#include "catch_amalgamated.hpp"

// System includes
// #include <iostream>

// Our includes
#include "CCZ4Geometry.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"

template <class data_t> struct vars_t
{
    data_t chi;
    Tensor<2, data_t> h;
    Tensor<1, data_t> Gamma;
};

TEST_CASE("CCZ4 Geometry")
{

    vars_t<double> vars;
    vars_t<Tensor<1, double>> d1;
    vars_t<Tensor<2, double>> d2;
    Tensor<1, double> Z_over_chi;

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
        CHECK_THAT(h_UU[i][j], Catch::Matchers::WithinAbs(h_UU_known[i][j],
                                                          test_threshold));
    }

    FOR (i, j, k)
    {
        INFO("chris.ULL[" << i << "][" << j << "][" << k << "]");
        CHECK_THAT(
            chris.ULL[i][j][k],
            Catch::Matchers::WithinAbs(chris_known[i][j][k], test_threshold));
    }

    FOR (i)
    {
        INFO("chris.contracted[" << i << "]");
        CHECK_THAT(chris.contracted[i],
                   Catch::Matchers::WithinAbs(chris_contracted_known[i],
                                              test_threshold));
    }

    FOR (i, j)
    {
        INFO("ricciZ.LL[" << i << "][" << j << "]");
        CHECK_THAT(ricciZ.LL[i][j], Catch::Matchers::WithinAbs(
                                        ricciZ_known[i][j], test_threshold));
    }

    CHECK_THAT(ricciZ.scalar,
               Catch::Matchers::WithinAbs(ricciZ_scalar_known, test_threshold));
}
