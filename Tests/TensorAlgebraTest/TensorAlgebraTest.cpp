/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// Doctest header
#include "doctest.h"

// Test header
#include "TensorAlgebraTest.hpp"

// Common includes
#include "doctestCLIArgs.hpp"

// AMReX includes
#include "AMReX.H"

// Our includes
#include "TensorAlgebra.hpp"

namespace
{
void check_tensor(const Tensor::Rank2 &tensor,
                  const Tensor::Rank2 &correct_tensor,
                  const std::string &test_name)
{
    FOR (i, j)
    {
        INFO(test_name << ": component [" << i << "][" << j << "]");
        CHECK(tensor(i, j) == correct_tensor(i, j));
    }
}

void check_tensor(const Tensor::Rank3 &tensor,
                  const Tensor::Rank3 &correct_tensor,
                  const std::string &test_name)
{
    FOR (i, j, k)
    {
        INFO(test_name << ": component [" << i << "][" << j << "][" << k
                       << "]");
        CHECK(tensor(i, j, k) == correct_tensor(i, j, k));
    }
}

void check_tensor(const Tensor::Rank4 &tensor,
                  const Tensor::Rank4 &correct_tensor,
                  const std::string &test_name)
{
    FOR (i, j, k, l)
    {
        INFO(test_name << ": component [" << i << "][" << j << "][" << k
                       << "][" << l << "]");
        CHECK(tensor(i, j, k, l) == correct_tensor(i, j, k, l));
    }
}
} // namespace

void run_tensor_algebra_test()
{
    int amrex_argc    = doctest::cli_args.argc();
    char **amrex_argv = doctest::cli_args.argv();
    // NOLINTNEXTLINE(bugprone-casting-through-void) // Open MPI triggers this
    amrex::Initialize(amrex_argc, amrex_argv);
    {
        using namespace TensorAlgebra;

        // Fill each symmetric Tensor with distinct values and set both
        // orderings of the symmetric indices in the expected general Tensor
        amrex::Real value = 0.;

        // Test make_general_from_symmetric for Sym12Rank2
        Tensor::Sym12Rank2 sym12_rank2{};
        Tensor::Rank2 correct_rank2{};
        FOR2_SYM(i, j)
        {
            value += 1.;
            sym12_rank2(i, j)   = value;
            correct_rank2(i, j) = value;
            correct_rank2(j, i) = value;
        }
        check_tensor(make_general_from_symmetric(sym12_rank2), correct_rank2,
                     "make_general_from_symmetric(Sym12Rank2)");

        // Test make_general_from_symmetric for Sym12Rank3
        Tensor::Sym12Rank3 sym12_rank3{};
        Tensor::Rank3 correct_rank3{};
        FOR (k)
        {
            FOR2_SYM(i, j)
            {
                value += 1.;
                sym12_rank3(i, j, k)   = value;
                correct_rank3(i, j, k) = value;
                correct_rank3(j, i, k) = value;
            }
        }
        check_tensor(make_general_from_symmetric(sym12_rank3), correct_rank3,
                     "make_general_from_symmetric(Sym12Rank3)");

        // Test make_general_from_symmetric for Sym23Rank3
        Tensor::Sym23Rank3 sym23_rank3{};
        FOR (i)
        {
            FOR2_SYM(j, k)
            {
                value += 1.;
                sym23_rank3(i, j, k)   = value;
                correct_rank3(i, j, k) = value;
                correct_rank3(i, k, j) = value;
            }
        }
        check_tensor(make_general_from_symmetric(sym23_rank3), correct_rank3,
                     "make_general_from_symmetric(Sym23Rank3)");

        // Test make_general_from_symmetric for Sym34Rank4
        Tensor::Sym34Rank4 sym34_rank4{};
        Tensor::Rank4 correct_rank4{};
        FOR (i, j)
        {
            FOR2_SYM(k, l)
            {
                value += 1.;
                sym34_rank4(i, j, k, l)   = value;
                correct_rank4(i, j, k, l) = value;
                correct_rank4(i, j, l, k) = value;
            }
        }
        check_tensor(make_general_from_symmetric(sym34_rank4), correct_rank4,
                     "make_general_from_symmetric(Sym34Rank4)");

        // Test make_general_from_symmetric for Sym12Sym34Rank4
        Tensor::Sym12Sym34Rank4 sym12_sym34_rank4{};
        FOR2_SYM(i, j)
        {
            FOR2_SYM(k, l)
            {
                value += 1.;
                sym12_sym34_rank4(i, j, k, l) = value;
                correct_rank4(i, j, k, l)     = value;
                correct_rank4(j, i, k, l)     = value;
                correct_rank4(i, j, l, k)     = value;
                correct_rank4(j, i, l, k)     = value;
            }
        }
        check_tensor(make_general_from_symmetric(sym12_sym34_rank4),
                     correct_rank4,
                     "make_general_from_symmetric(Sym12Sym34Rank4)");
    }
    amrex::Finalize();
}
