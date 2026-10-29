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

using namespace amrex::literals;

namespace
{
constexpr int ulp = 15; // units in the last place
constexpr amrex::Real real_epsilon =
    std::numeric_limits<amrex::Real>::epsilon();

void check_tensor(const Tensor::Rank2 &tensor,
                  const Tensor::Rank2 &correct_tensor,
                  const std::string &test_name)
{
    FOR (i, j)
    {
        INFO(test_name << ": component [" << i << "][" << j << "]");
        CHECK(
            tensor(i, j) ==
            doctest::Approx(correct_tensor(i, j)).epsilon(ulp * real_epsilon));
    }
}

void check_tensor(const Tensor::Sym12Rank2 &tensor,
                  const Tensor::Sym12Rank2 &correct_tensor,
                  const std::string &test_name)
{
    FOR (i, j)
    {
        INFO(test_name << ": component [" << i << "][" << j << "]");
        CHECK(
            tensor(i, j) ==
            doctest::Approx(correct_tensor(i, j)).epsilon(ulp * real_epsilon));
    }
}

void check_vector(const Tensor::Rank1 &vector,
                  const Tensor::Rank1 &correct_vector,
                  const std::string &test_name)
{
    FOR (i)
    {
        INFO(test_name << ": component [" << i << "]");
        CHECK(vector(i) ==
              doctest::Approx(correct_vector(i)).epsilon(ulp * real_epsilon));
    }
}

// Check that rotation is a proper rotation (R^T R = I and det R = 1) which
// maps the direction of origin onto the direction of destination
// NOLINTBEGIN(bugprone-easily-swappable-parameters)
void check_rotation(const Tensor::Rank2 &rotation, const Tensor::Rank1 &origin,
                    const Tensor::Rank1 &destination,
                    const std::string &test_name)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    amrex::Real norm_origin = std::sqrt(
        origin(0) * origin(0) + origin(1) * origin(1) + origin(2) * origin(2));
    amrex::Real norm_dest = std::sqrt(destination(0) * destination(0) +
                                      destination(1) * destination(1) +
                                      destination(2) * destination(2));
    Tensor::Rank1 origin_hat{};
    Tensor::Rank1 dest_hat{};
    FOR (i)
    {
        origin_hat(i) = origin(i) / norm_origin;
        dest_hat(i)   = destination(i) / norm_dest;
    }
    check_vector(
        CoordinateTransformations::transform_vector(origin_hat, rotation),
        dest_hat, test_name + " (R origin = destination)");

    Tensor::Rank2 identity{};
    Tensor::Rank2 rotation_T_rotation{};
    FOR (i, j)
    {
        identity(i, j)            = (i == j) ? 1. : 0.;
        rotation_T_rotation(i, j) = 0.;
        FOR (k)
        {
            rotation_T_rotation(i, j) += rotation(k, i) * rotation(k, j);
        }
    }
    check_tensor(rotation_T_rotation, identity, test_name + " (R^T R = I)");

    INFO(test_name << ": determinant");
    CHECK(TensorAlgebra::compute_determinant(rotation) ==
          doctest::Approx(1.).epsilon(ulp * real_epsilon));
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
        amrex::Real rho2        = std::max(x * x + y * y, 1e-12_rt);
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
        Tensor::Rank2 jac           = spherical_jacobian(x, y, z);
        Tensor::Rank2 inv_jac       = inverse_spherical_jacobian(x, y, z);
        Tensor::Rank2 inv_jac_check = compute_inverse(jac);
        check_tensor(inv_jac, inv_jac_check, "inverse_jacobian");

        // Test tensor transformations
        Tensor::Sym12Rank2 Mij_cart{0.};
        Mij_cart(0, 0) = 1.;
        Mij_cart(1, 1) = 1.;
        Mij_cart(2, 2) = 1.;

        Tensor::Sym12Rank2 Mij_spher{0.};
        Mij_spher(0, 0) = 1.;
        Mij_spher(1, 1) = r * r;
        Mij_spher(2, 2) = r2sin2theta;

        // Test cartesian_to_spherical_LL
        Tensor::Sym12Rank2 Mij_spher_check{};
        Mij_spher_check = cartesian_to_spherical_LL(Mij_cart, x, y, z);
        check_tensor(Mij_spher_check, Mij_spher, "cartesian_to_spherical_LL");

        // Test spherical_to_cartesian_LL
        Tensor::Sym12Rank2 Mij_cart_check{};
        Mij_cart_check = spherical_to_cartesian_LL(Mij_spher, x, y, z);
        check_tensor(Mij_cart_check, Mij_cart, "spherical_to_cartesian_LL");

        // Test cartesian_to_spherical_UU
        Tensor::Sym12Rank2 Mij_spher_UU{};
        Tensor::Sym12Rank2 Mij_spher_UU_check{};
        Mij_spher_UU_check =
            cartesian_to_spherical_UU(compute_inverse_sym(Mij_cart), x, y, z);
        Mij_spher_UU = compute_inverse_sym(Mij_spher);
        check_tensor(Mij_spher_UU_check, Mij_spher_UU,
                     "cartesian_to_spherical_UU");

        // Test spherical_to_cartesian_UU
        Tensor::Sym12Rank2 Mij_cart_UU{};
        Tensor::Sym12Rank2 Mij_cart_UU_check{};
        Mij_cart_UU_check =
            spherical_to_cartesian_UU(compute_inverse_sym(Mij_spher), x, y, z);
        Mij_cart_UU = compute_inverse_sym(Mij_cart);
        check_tensor(Mij_cart_UU_check, Mij_cart_UU,
                     "spherical_to_cartesian_UU");

        // Test the symmetric transformations on a generic symmetric tensor
        // with non-zero off-diagonal components (the flat metric above is
        // diagonal in both coordinate systems so cannot detect errors in the
        // off-diagonal components). Compare against the full Rank2
        // transformations which compute the same contraction.
        Tensor::Sym12Rank2 Nij{};
        Nij(0, 0) = 1.3;
        Nij(0, 1) = 0.7;
        Nij(0, 2) = -0.4;
        Nij(1, 1) = 2.1;
        Nij(1, 2) = 0.9;
        Nij(2, 2) = 1.6;

        const Tensor::Rank2 Nij_full = make_general_from_symmetric(Nij);

        check_tensor(make_general_from_symmetric(
                         spherical_to_cartesian_LL(Nij, x, y, z)),
                     transform_tensor_LL(Nij_full, jac),
                     "spherical_to_cartesian_LL (generic)");
        check_tensor(make_general_from_symmetric(
                         spherical_to_cartesian_UU(Nij, x, y, z)),
                     transform_tensor_UU(Nij_full, inv_jac),
                     "spherical_to_cartesian_UU (generic)");
        check_tensor(make_general_from_symmetric(
                         cartesian_to_spherical_LL(Nij, x, y, z)),
                     transform_tensor_LL(Nij_full, inv_jac),
                     "cartesian_to_spherical_LL (generic)");
        check_tensor(make_general_from_symmetric(
                         cartesian_to_spherical_UU(Nij, x, y, z)),
                     transform_tensor_UU(Nij_full, jac),
                     "cartesian_to_spherical_UU (generic)");

        // Test vector transformations
        Tensor::Rank1 si_cart{};
        si_cart(0) = x / r;
        si_cart(1) = y / r;
        si_cart(2) = z / r;

        Tensor::Rank1 si_spher{};
        si_spher(0) = 1.0;
        si_spher(1) = 0.0;
        si_spher(2) = 0.0;

        // Test cartesian_to_spherical_U
        Tensor::Rank1 si_spher_U_check{};
        si_spher_U_check = cartesian_to_spherical_U(si_cart, x, y, z);
        check_vector(si_spher_U_check, si_spher, "cartesian_to_spherical_U");

        // Test spherical_to_cartesian_U
        Tensor::Rank1 si_cart_U_check{};
        si_cart_U_check = spherical_to_cartesian_U(si_spher, x, y, z);
        check_vector(si_cart_U_check, si_cart, "spherical_to_cartesian_U");

        // Test cartesian_to_spherical_L
        Tensor::Rank1 si_spher_L_check{};
        si_spher_L_check = cartesian_to_spherical_L(si_cart, x, y, z);
        check_vector(si_spher_L_check, si_spher, "cartesian_to_spherical_L");

        // Test spherical_to_cartesian_L
        Tensor::Rank1 si_cart_L_check{};
        si_cart_L_check = spherical_to_cartesian_L(si_spher, x, y, z);
        check_vector(si_cart_L_check, si_cart, "spherical_to_cartesian_L");

        // Test area_element_sphere
        amrex::Real area_element       = r * sqrt(rho2);
        amrex::Real area_element_check = area_element_sphere(Mij_spher);
        CHECK(area_element ==
              doctest::Approx(area_element_check).epsilon(ulp * real_epsilon));

        // Test rotation_matrix(axis, cos_angle)
        Tensor::Rank2 identity{};
        FOR (i, j)
        {
            identity(i, j) = (i == j) ? 1. : 0.;
        }

        // rotation by 90 degrees about z
        Tensor::Rank2 rotation_z_90{};
        rotation_z_90(0, 1) = -1.;
        rotation_z_90(1, 0) = 1.;
        rotation_z_90(2, 2) = 1.;
        check_tensor(rotation_matrix({0., 0., 1.}, 0.), rotation_z_90,
                     "rotation_matrix (90 degrees about z)");

        // rotation by 0 degrees
        check_tensor(rotation_matrix({0., 0., 1.}, 1.), identity,
                     "rotation_matrix (0 degrees)");

        // rotation about a generic axis leaves the axis unchanged
        const amrex::Real norm_123 = std::sqrt(14.);
        Tensor::Rank1 axis_123 = {1. / norm_123, 2. / norm_123, 3. / norm_123};
        check_rotation(rotation_matrix(axis_123, 0.3), axis_123, axis_123,
                       "rotation_matrix (generic axis)");

        // Test rotation_matrix(origin, destination)
        auto check_rotation_between = [](const Tensor::Rank1 &origin,
                                         const Tensor::Rank1 &destination,
                                         const std::string &test_name)
        {
            check_rotation(rotation_matrix(origin, destination), origin,
                           destination, test_name);
        };

        // KerrBH default spin direction onto the z axis
        check_rotation_between({0.1, 0.1, 0.}, {0., 0., 1.},
                               "rotation_matrix (KerrBH spin direction)");
        // generic vectors which are not normalized
        check_rotation_between({1., 2., 3.}, {-2., 0.5, 1.},
                               "rotation_matrix (generic)");

        // parallel vectors give the identity
        check_tensor(rotation_matrix({1., 2., 3.}, {2., 4., 6.}), identity,
                     "rotation_matrix (parallel)");

        // anti-parallel vectors give a rotation by 180 degrees
        check_rotation_between({0., 0., -1.}, {0., 0., 1.},
                               "rotation_matrix (anti-parallel z)");
        check_rotation_between({1., 2., 3.}, {-1., -2., -3.},
                               "rotation_matrix (anti-parallel generic)");
        check_rotation_between({1., 0., 0.}, {-1., 0., 0.},
                               "rotation_matrix (anti-parallel x)");
    }
    amrex::Finalize();
}
