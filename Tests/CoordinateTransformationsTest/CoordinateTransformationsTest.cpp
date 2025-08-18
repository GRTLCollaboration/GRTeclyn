/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// Doctest header
#include "doctest.h"

// Test include
#include "CoordinateTransformationsTest.hpp"

// Common includes
#include "doctestCLIArgs.hpp"

// AMReX includes
#include "AMReX.H"
#include "AMReX_IntVect.H"

// Other includes
#include <limits>

// Our includes
#include "CoordinateTransformations.hpp"
#include "Coordinates.hpp"
#include "TensorAlgebra.hpp"

namespace
{
constexpr int ulp             = 15; // units in the last place
constexpr double real_epsilon = std::numeric_limits<amrex::Real>::epsilon();

void check_tensor(const Tensor<2, amrex::Real> &tensor,
                  const Tensor<2, amrex::Real> &correct_tensor,
                  const std::string &test_name)
{
    FOR (i, j)
    {
        INFO(test_name << ": component [" << i << "][" << j << "]");
        CHECK(
            tensor[i][j] ==
            doctest::Approx(correct_tensor[i][j]).epsilon(ulp * real_epsilon));
    }
}

void check_vector(const Tensor<1, amrex::Real> &vector,
                  const Tensor<1, amrex::Real> &correct_vector,
                  const std::string &test_name)
{
    FOR (i)
    {
        INFO(test_name << ": component [" << i << "]");
        CHECK(vector[i] ==
              doctest::Approx(correct_vector[i]).epsilon(ulp * real_epsilon));
    }
}
} // namespace

void run_coordinate_transformations_test()
{
    int amrex_argc    = doctest::cli_args.argc();
    char **amrex_argv = doctest::cli_args.argv();
    // NOLINTNEXTLINE(bugprone-casting-through-void) // Open MPI triggers this
    amrex::Initialize(amrex_argc, amrex_argv);
    {
        const amrex::Real dx = 0.1;
        amrex::IntVect iv{1, 2, 3};

        Coordinates coords(iv, dx);
        const amrex::Real x     = coords.x;
        const amrex::Real y     = coords.y;
        const amrex::Real z     = coords.z;
        const amrex::Real r     = coords.get_radius();
        amrex::Real rho2        = std::max(x * x + y * y, 1e-12);
        amrex::Real r2sin2theta = rho2;

        /* for debugging
        std::cout << "x " << x << std::endl;
        std::cout << "y " << y << std::endl;
        std::cout << "z " << z << std::endl;
        std::cout << "r " << r << std::endl;
        */

        using namespace TensorAlgebra;
        using namespace CoordinateTransformations;

        // Test if inv_jac is really the inverse of the jacobian
        Tensor<2, amrex::Real> jac     = spherical_jacobian(x, y, z);
        Tensor<2, amrex::Real> inv_jac = inverse_spherical_jacobian(x, y, z);
        Tensor<2, amrex::Real> inv_jac_check = compute_inverse(jac);
        check_tensor(inv_jac, inv_jac_check, "inverse_jacobian");

        // Test tensor transformations
        Tensor<2, amrex::Real> Mij_cart;
        FOR (i, j)
        {
            Mij_cart[i][j] = 0.;
        }
        Mij_cart[0][0] = 1.;
        Mij_cart[1][1] = 1.;
        Mij_cart[2][2] = 1.;

        Tensor<2, amrex::Real> Mij_spher;
        FOR (i, j)
        {
            Mij_spher[i][j] = 0.;
        }
        Mij_spher[0][0] = 1.;
        Mij_spher[1][1] = r * r;
        Mij_spher[2][2] = r2sin2theta;

        // Test cartesian_to_spherical_LL
        Tensor<2, amrex::Real> Mij_spher_check;
        Mij_spher_check = cartesian_to_spherical_LL(Mij_cart, x, y, z);
        check_tensor(Mij_spher_check, Mij_spher, "cartesian_to_spherical_LL");

        // Test spherical_to_cartesian_LL
        Tensor<2, amrex::Real> Mij_cart_check;
        Mij_cart_check = spherical_to_cartesian_LL(Mij_spher, x, y, z);
        check_tensor(Mij_cart_check, Mij_cart, "spherical_to_cartesian_LL");

        // Test cartesian_to_spherical_UU
        Tensor<2, amrex::Real> Mij_spher_UU;
        Tensor<2, amrex::Real> Mij_spher_UU_check;
        Mij_spher_UU_check =
            cartesian_to_spherical_UU(compute_inverse_sym(Mij_cart), x, y, z);
        Mij_spher_UU = compute_inverse_sym(Mij_spher);
        check_tensor(Mij_spher_UU_check, Mij_spher_UU,
                     "cartesian_to_spherical_UU");

        // Test spherical_to_cartesian_UU
        Tensor<2, amrex::Real> Mij_cart_UU;
        Tensor<2, amrex::Real> Mij_cart_UU_check;
        Mij_cart_UU_check =
            spherical_to_cartesian_UU(compute_inverse_sym(Mij_spher), x, y, z);
        Mij_cart_UU = compute_inverse_sym(Mij_cart);
        check_tensor(Mij_cart_UU_check, Mij_cart_UU,
                     "spherical_to_cartesian_UU");

        // Test vector transformations
        Tensor<1, amrex::Real> si_cart;
        si_cart[0] = x / r;
        si_cart[1] = y / r;
        si_cart[2] = z / r;

        Tensor<1, amrex::Real> si_spher;
        si_spher[0] = 1.0;
        si_spher[1] = 0.0;
        si_spher[2] = 0.0;

        // Test cartesian_to_spherical_U
        Tensor<1, amrex::Real> si_spher_U_check;
        si_spher_U_check = cartesian_to_spherical_U(si_cart, x, y, z);
        check_vector(si_spher_U_check, si_spher, "cartesian_to_spherical_U");

        // Test spherical_to_cartesian_U
        Tensor<1, amrex::Real> si_cart_U_check;
        si_cart_U_check = spherical_to_cartesian_U(si_spher, x, y, z);
        check_vector(si_cart_U_check, si_cart, "spherical_to_cartesian_U");

        // Test cartesian_to_spherical_L
        Tensor<1, amrex::Real> si_spher_L_check;
        si_spher_L_check = cartesian_to_spherical_L(si_cart, x, y, z);
        check_vector(si_spher_L_check, si_spher, "cartesian_to_spherical_L");

        // Test spherical_to_cartesian_L
        Tensor<1, amrex::Real> si_cart_L_check;
        si_cart_L_check = spherical_to_cartesian_L(si_spher, x, y, z);
        check_vector(si_cart_L_check, si_cart, "spherical_to_cartesian_L");

        // Test area_element_sphere
        amrex::Real area_element       = r * sqrt(rho2);
        amrex::Real area_element_check = area_element_sphere(Mij_spher);
        CHECK(area_element ==
              doctest::Approx(area_element_check).epsilon(ulp * real_epsilon));
    }
    amrex::Finalize();
}
