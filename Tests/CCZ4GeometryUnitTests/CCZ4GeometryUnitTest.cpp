/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// Doctest header
#include "doctest.h"
#include "doctestCLIArgs.hpp"

// Test header
#include "CCZ4GeometryUnitTest.hpp"

// AMReX headers
#include "AMReX.H"
#include "AMReX_FArrayBox.H"

// System includes
#include <iostream>

// Our includes
#include "CCZ4D1Vars.hpp"
#include "CCZ4D2Vars.hpp"
#include "CCZ4Geometry.hpp"
#include "CCZ4Vars.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
void run_ccz4_geometry_unit_tests()
{

    int amrex_argc    = doctest::cli_args.argc();
    char **amrex_argv = doctest::cli_args.argv();
    // NOLINTNEXTLINE(bugprone-casting-through-void) // Open MPI triggers this
    amrex::Initialize(amrex_argc, amrex_argv);
    {
        amrex::Real chi = 0.0;
        Tensor<1, amrex::Real> Gamma;
        Tensor<2, amrex::Real> h;
        CCZ4D1Vars d1;
        CCZ4D2Vars d2;
        Tensor<1, amrex::Real> Z_over_chi;

// Including the auto generated file with values
#include "CCZ4GeometryMathematicaValues.hpp"

        // Get a CCZ4Vars object and store the test values
        const amrex::IntVect iv_zeros(0, 0, 0);
        const amrex::Box box(iv_zeros, iv_zeros);
        amrex::FArrayBox in_fab{box, NUM_CCZ4_VARS, amrex::The_Managed_Arena()};
        in_fab.setVal(0.0);
        const amrex::Array4<amrex::Real> &in_array = in_fab.array();
        in_array(0, 0, 0, c_chi)                   = chi;
        FOR (i)
        {
            in_array(iv_zeros, c_Gamma1 + i) = Gamma[i];
            FOR (j)
            {
                in_array(iv_zeros, var_idx(c_h11, i, j)) = h[i][j];
            }
        }

        const amrex::CellData<const amrex::Real> &cell_data =
            in_array.cellData(0, 0, 0);
        ConstCCZ4Vars vars(cell_data);

        auto h_UU  = CCZ4Geometry::compute_inverse_metric(vars);
        auto chris = TensorAlgebra::compute_christoffel(d1.h, h_UU);

        auto ricciZ = CCZ4Geometry::compute_ricci_Z(vars, d1, d2.chi, d2.h,
                                                    h_UU, chris, Z_over_chi);

        double test_threshold = 1e-14;

        // Compare
        FOR (i, j)
        {
            INFO("h_UU[" << i << "][" << j << "]");
            CHECK(h_UU[i][j] ==
                  doctest::Approx(h_UU_known[i][j]).epsilon(test_threshold));
        }

        FOR (i, j, k)
        {
            INFO("chris.ULL[" << i << "][" << j << "][" << k << "]");
            CHECK(
                chris.ULL[i][j][k] ==
                doctest::Approx(chris_known[i][j][k]).epsilon(test_threshold));
        }

        FOR (i)
        {
            INFO("chris.contracted[" << i << "]");
            CHECK(chris.contracted[i] ==
                  doctest::Approx(chris_contracted_known[i])
                      .epsilon(test_threshold));
        }

        FOR (i, j)
        {
            INFO("ricciZ.LL[" << i << "][" << j << "]");
            CHECK(ricciZ.LL[i][j] ==
                  doctest::Approx(ricciZ_known[i][j]).epsilon(test_threshold));
        }

        CHECK(ricciZ.scalar ==
              doctest::Approx(ricciZ_scalar_known).epsilon(test_threshold));
    }
    amrex::Finalize();
}
