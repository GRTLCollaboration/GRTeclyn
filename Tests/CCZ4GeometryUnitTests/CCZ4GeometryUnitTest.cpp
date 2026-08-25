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
#include "CCZ4Geometry.hpp"
#include "CCZ4Vars.hpp"
#include "DimensionDefinitions.hpp"

// NOLINTBEGIN(readability-function-cognitive-complexity)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
compute_ccz4_test_geometry(const amrex::Array4<amrex::Real> &a_array,
                           const amrex::IntVect &a_iv,
                           const amrex::Array4<amrex::Real> &a_geometry_array)
{
    amrex::Real chi = 0.0;
    Tensor::Rank1 Gamma{};
    Tensor::Rank2 h{};
    Tensor::Rank1 Z_over_chi{};
    Tensor::Rank1 d1_chi{};
    Tensor::Rank2 d1_Gamma{};
    Tensor::Sym12Rank3 d1_h{};
    Tensor::Sym12Sym34Rank4 d2_h{};
    Tensor::Sym12Rank2 d2_chi{};
// Including the auto generated file with values
// Run it with "python CCZ4GeometryGenerateExpectedValues.py"
#define CCZ4_GEOMETRY_INPUT_VALUES
#include "CCZ4GeometryExpectedValues.hpp"
#undef CCZ4_GEOMETRY_INPUT_VALUES

    a_array(a_iv, c_chi) = chi;
    FOR (i)
    {

        a_array(a_iv, c_Gamma1 + i) = Gamma(i);
        FOR (j)
        {
            a_array(a_iv, sym_var_idx(c_h11, i, j)) = h(i, j);
        }
    }

    const amrex::CellData<const amrex::Real> &cell_data =
        a_array.cellData(a_iv[0], a_iv[1], a_iv[2]);
    CCZ4Vars vars(cell_data);

    auto h_UU  = CCZ4Geometry::compute_inverse_metric(vars);
    auto chris = CCZ4Geometry::compute_christoffel(d1_h, h_UU);
    auto phys_chris =
        CCZ4Geometry::compute_phys_chris(vars, d1_chi, h_UU, chris.ULL);
    auto ricciZ = CCZ4Geometry::compute_ricci_Z(
        vars, d1_chi, d1_Gamma, d1_h, d2_h, d2_chi, h_UU, chris, Z_over_chi);
    const amrex::Real dZ_coeff = 1.0;
    auto ricciZ_general        = CCZ4Geometry::compute_ricci_Z_general(
        vars, d1_chi, d1_Gamma, d1_h, d2_chi, d2_h, h_UU, chris, dZ_coeff);

    int vars_counter = 0;
    FOR (i, j)
    {
        a_geometry_array(a_iv, vars_counter) = h_UU(i, j);
        ++vars_counter;
    }
    FOR (i, j, k)
    {

        a_geometry_array(a_iv, vars_counter) = chris.ULL(i, j, k);
        ++vars_counter;
    }

    FOR (i)
    {

        a_geometry_array(a_iv, vars_counter) = chris.contracted(i);
        ++vars_counter;
    }

    FOR (i, j, k)
    {

        a_geometry_array(a_iv, vars_counter) = phys_chris(i, j, k);
        ++vars_counter;
    }

    FOR (i, j)
    {
        a_geometry_array(a_iv, vars_counter) = ricciZ.LL(i, j);
        ++vars_counter;
    }

    a_geometry_array(a_iv, vars_counter++) = ricciZ.scalar;

    FOR (i, j)
    {
        a_geometry_array(a_iv, vars_counter) = ricciZ_general.LL(i, j);
        ++vars_counter;
    }

    a_geometry_array(a_iv, vars_counter) = ricciZ_general.scalar;
}

void run_ccz4_geometry_unit_tests()
{
#define CCZ4_GEOMETRY_EXPECTED_VALUES
#include "CCZ4GeometryExpectedValues.hpp"
#undef CCZ4_GEOMETRY_EXPECTED_VALUES

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

        amrex::Real test_threshold = 1e-14;

        const amrex::CellData<const amrex::Real> &geometry_test_cell_data =
            geometry_array.cellData(0, 0, 0);

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

        FOR (i, j, k)
        {
            INFO("phys_chris.ULL[" << i << "][" << j << "][" << k << "]");
            CHECK(geometry_array(iv_zeros, vars_counter) ==
                  doctest::Approx(chris_phys_known[i][j][k])
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

        FOR (i, j)
        {
            INFO("ricciZ_general.LL[" << i << "][" << j << "]");
            CHECK(geometry_array(iv_zeros, vars_counter) ==
                  doctest::Approx(ricciZ_general_known[i][j])
                      .epsilon(test_threshold));
            ++vars_counter;
        }

        CHECK(geometry_array(iv_zeros, vars_counter) ==
              doctest::Approx(ricciZ_general_scalar_known)
                  .epsilon(test_threshold));
        ++vars_counter;
    }
    amrex::Finalize();
}

// NOLINTEND(readability-function-cognitive-complexity)
