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
#include "AMReX_RealVect.H"

// System includes
#include <iostream>

// Our includes
#include "CCZ4D1Vars.hpp"
#include "CCZ4D2Vars.hpp"
#include "CCZ4Geometry.hpp"
#include "CCZ4Vars.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
compute_ccz4_test_geometry(const amrex::Array4<amrex::Real> &a_array,
                           const amrex::IntVect &a_iv,
                           const amrex::Array4<amrex::Real> &a_geometry_array)
{
    amrex::Real chi = 0.0;
    Tensor<1, amrex::Real> Gamma;
    Tensor<2, amrex::Real> h;
    CCZ4D1Vars d1;
    CCZ4D2Vars d2;
    Tensor<1, amrex::Real> Z_over_chi;

// Including the auto generated file with values
#include "CCZ4GeometryMathematicaValues.hpp"

    a_array(a_iv, c_chi) = chi;
    FOR (i)
    {
        a_array(a_iv, c_Gamma1 + i) = Gamma[i];
        FOR (j)
        {
            a_array(a_iv, VAR_IDX(c_h11, i, j)) = h[i][j];
        }
    }

    const amrex::CellData<const amrex::Real> &cell_data =
        a_array.cellData(a_iv[0], a_iv[1], a_iv[2]);
    CCZ4Vars vars(cell_data);

    auto h_UU   = CCZ4Geometry::compute_inverse_metric(vars);
    auto chris  = CCZ4Geometry::compute_christoffel(d1, h_UU);
    auto ricciZ = CCZ4Geometry::compute_ricci_Z(vars, d1, d2.chi, d2.h, h_UU,
                                                chris, Z_over_chi);

    int vars_counter = 0;
    FOR (i, j)
    {
        a_geometry_array(a_iv, vars_counter) = h_UU[i][j];
        ++vars_counter;
    }
    FOR (i, j, k)
    {
        a_geometry_array(a_iv, vars_counter) = chris.ULL[i][j][k];
        ++vars_counter;
    }

    FOR (i)
    {
        a_geometry_array(a_iv, vars_counter) = chris.contracted[i];
        ++vars_counter;
    }

    FOR (i, j)
    {
        a_geometry_array(a_iv, vars_counter) = ricciZ.LL[i][j];
        ++vars_counter;
    }

    a_geometry_array(a_iv, vars_counter) = ricciZ.scalar;
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
void run_ccz4_geometry_unit_tests()
{

    int amrex_argc    = doctest::cli_args.argc();
    char **amrex_argv = doctest::cli_args.argv();
    // NOLINTNEXTLINE(bugprone-casting-through-void) // Open MPI triggers this
    amrex::Initialize(amrex_argc, amrex_argv);
    {
        // Get a CCZ4Vars object and store the test values
        const amrex::IntVect iv_zeros(0, 0, 0);
        const amrex::Box box(iv_zeros, iv_zeros);
        amrex::FArrayBox in_fab{box, NUM_CCZ4_VARS, amrex::The_Managed_Arena()};
        amrex::FArrayBox geometry_fab{box, NUM_GEOMETRY_TEST_VARS,
                                      amrex::The_Managed_Arena()};
        const amrex::Array4<amrex::Real> &in_array = in_fab.array();
        const auto &geometry_array                 = geometry_fab.array();

        amrex::ParallelFor(box,
                           [=] AMREX_GPU_DEVICE(int ix, int iy, int iz)
                           {
                               const amrex::IntVect iv{ix, iy, iz};
                               compute_ccz4_test_geometry(in_array, iv,
                                                          geometry_array);
                           });

        amrex::Gpu::streamSynchronize();

        double test_threshold = 1e-14;

        const amrex::CellData<const amrex::Real> &geometry_test_cell_data =
            geometry_array.cellData(0, 0, 0);

// Including the auto generated file with expected values
#include "CCZ4GeometryMathematicaExpectedValues.hpp"

        int vars_counter = 0;
        // Compare
        FOR (i, j)
        {
            INFO("h_UU[" << i << "][" << j << "]");
            CHECK(geometry_array(iv_zeros, vars_counter) ==
                  doctest::Approx(h_UU_known[i][j]).epsilon(test_threshold));
            ++vars_counter;
        }

        FOR (i, j, k)
        {
            INFO("chris.ULL[" << i << "][" << j << "][" << k << "]");
            CHECK(
                geometry_array(iv_zeros, vars_counter) ==
                doctest::Approx(chris_known[i][j][k]).epsilon(test_threshold));
            ++vars_counter;
        }

        FOR (i)
        {
            INFO("chris.contracted[" << i << "]");
            CHECK(geometry_array(iv_zeros, vars_counter) ==
                  doctest::Approx(chris_contracted_known[i])
                      .epsilon(test_threshold));
            ++vars_counter;
        }

        FOR (i, j)
        {
            INFO("ricciZ.LL[" << i << "][" << j << "]");
            CHECK(geometry_array(iv_zeros, vars_counter) ==
                  doctest::Approx(ricciZ_known[i][j]).epsilon(test_threshold));
            ++vars_counter;
        }

        CHECK(geometry_array(iv_zeros, vars_counter) ==
              doctest::Approx(ricciZ_scalar_known).epsilon(test_threshold));
        ++vars_counter;
    }
    amrex::Finalize();
}
