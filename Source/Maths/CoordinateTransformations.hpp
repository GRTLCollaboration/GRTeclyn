/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef COORDINATETRANSFORMATIONS_HPP_
#define COORDINATETRANSFORMATIONS_HPP_

#include "DimensionDefinitions.hpp"
#include "TensorAlgebra.hpp"

#include <cmath>

using namespace amrex::literals;

namespace CoordinateTransformations
{

// Jacobian transformation matrix

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
static Tensor::Rank2 spherical_jacobian(const amrex::Real x,
                                        const amrex::Real y,
                                        const amrex::Real z)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    // calculate useful position quantities
    amrex::Real rho2 = x * x + y * y;
    rho2             = std::max(rho2, 1e-12_rt);
    amrex::Real rho  = std::sqrt(rho2);
    amrex::Real r2   = x * x + y * y + z * z;
    r2               = std::max(r2, 1e-12_rt);
    amrex::Real r    = std::sqrt(r2);

    // And the sines and cosines of phi and theta
    amrex::Real cos_phi = x / rho;
    amrex::Real sin_phi = y / rho;

    // derivatives for jacobian matrix - drdx etc
    Tensor::Rank2 jac{};
    jac(0, 0) = x / r;
    jac(1, 0) = cos_phi * z / r2;
    jac(2, 0) = -y / rho2;
    jac(0, 1) = y / r;
    jac(1, 1) = sin_phi * z / r2;
    jac(2, 1) = x / rho2;
    jac(0, 2) = z / r;
    jac(1, 2) = -rho / r2;
    jac(2, 2) = 0.0;
    return jac;
}

// Inverse Jacobian

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
static Tensor::Rank2 inverse_spherical_jacobian(const amrex::Real x,
                                                const amrex::Real y,
                                                const amrex::Real z)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    // calculate useful position quantities
    amrex::Real rho2 = x * x + y * y;
    amrex::Real rho  = std::sqrt(rho2);
    rho              = std::max(rho, 1e-6_rt);
    amrex::Real r2   = x * x + y * y + z * z;
    amrex::Real r    = std::sqrt(r2);
    r                = std::max(r, 1e-6_rt);

    // And the sines and cosines of phi and theta
    amrex::Real cos_phi = x / rho;
    amrex::Real sin_phi = y / rho;

    // derivatives for inverse jacobian matrix - drdx etc
    Tensor::Rank2 inv_jac{};
    inv_jac(0, 0) = x / r;
    inv_jac(1, 0) = y / r;
    inv_jac(2, 0) = z / r;
    inv_jac(0, 1) = z * cos_phi;
    inv_jac(1, 1) = z * sin_phi;
    inv_jac(2, 1) = -rho;
    inv_jac(0, 2) = -y;
    inv_jac(1, 2) = x;
    inv_jac(2, 2) = 0.0;
    return inv_jac;
}

// Convert a Tensor (with two lower indices) in spherical coords to cartesian
// coords

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
static Tensor::Sym12Rank2
spherical_to_cartesian_LL(const Tensor::Sym12Rank2 &spherical_g,
                          const amrex::Real x, const amrex::Real y,
                          const amrex::Real z)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    Tensor::Sym12Rank2 cartesian_g{0.};

    // derivatives for jacobian matrix - drdx etc
    Tensor::Rank2 jac = spherical_jacobian(x, y, z);

    // Convert the Tensor to cartesian coords
    FOR (i, j)
    {
        FOR (k, l)
        {
            cartesian_g(i, j) += spherical_g(k, l) * jac(k, i) * jac(l, j);
        }
    }
    return cartesian_g;
}

// Convert a Tensor (with two upper indices) in spherical coords to cartesian
// coords

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
static Tensor::Sym12Rank2
spherical_to_cartesian_UU(const Tensor::Sym12Rank2 &spherical_g_UU,
                          const amrex::Real x, const amrex::Real y,
                          const amrex::Real z)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    Tensor::Sym12Rank2 cartesian_g_UU{0.};

    // derivatives for jacobian matrix - drdx etc
    Tensor::Rank2 inv_jac = inverse_spherical_jacobian(x, y, z);

    // Convert the Tensor to cartesian coords
    FOR (i, j)
    {
        FOR (k, l)
        {
            cartesian_g_UU(i, j) +=
                spherical_g_UU(k, l) * inv_jac(i, k) * inv_jac(j, l);
        }
    }
    return cartesian_g_UU;
}

// Convert a Tensor (with two lower indices) in cartesian coords to spherical
// coords

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
static Tensor::Sym12Rank2
cartesian_to_spherical_LL(const Tensor::Sym12Rank2 &cartesian_g,
                          const amrex::Real x, const amrex::Real y,
                          const amrex::Real z)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    Tensor::Sym12Rank2 spherical_g{0.};

    // derivatives for inverse jacobian matrix - drdx etc
    Tensor::Rank2 inv_jac = inverse_spherical_jacobian(x, y, z);

    // Convert the Tensor to spherical coords
    FOR (i, j)
    {
        FOR (k, l)
        {
            spherical_g(i, j) +=
                cartesian_g(k, l) * inv_jac(k, i) * inv_jac(l, j);
        }
    }
    return spherical_g;
}

// Convert a Tensor (with two upper indices) in cartesian coords to spherical
// coords

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
static Tensor::Sym12Rank2
cartesian_to_spherical_UU(const Tensor::Sym12Rank2 &cartesian_g_UU,
                          amrex::Real x, amrex::Real y, amrex::Real z)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    Tensor::Sym12Rank2 spherical_g_UU{0.};

    // derivatives for jacobian matrix - drdx etc
    Tensor::Rank2 jac = spherical_jacobian(x, y, z);

    // Convert the Tensor to spherical coords
    FOR (i, j)
    {
        FOR (k, l)
        {
            spherical_g_UU(i, j) +=
                cartesian_g_UU(k, l) * jac(i, k) * jac(j, l);
        }
    }
    return spherical_g_UU;
}

// Convert a vector (with one upper index) in spherical coords to cartesian
// coords

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
Tensor::Rank1 spherical_to_cartesian_U(const Tensor::Rank1 &spherical_v_U,
                                       amrex::Real x, amrex::Real y,
                                       amrex::Real z)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    Tensor::Rank1 cartesian_v_U{0., 0., 0.};

    // derivatives for inverse jacobian matrix - drdx etc
    Tensor::Rank2 inv_jac = inverse_spherical_jacobian(x, y, z);

    // transform the vector to cartesian coords
    FOR (i)
    {
        FOR (j)
        {
            cartesian_v_U(i) += inv_jac(i, j) * spherical_v_U(j);
        }
    }
    return cartesian_v_U;
}

// Convert a vector (with one lower index) in spherical coords to cartesian
// coords

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
Tensor::Rank1 spherical_to_cartesian_L(const Tensor::Rank1 &spherical_v_L,
                                       amrex::Real x, amrex::Real y,
                                       amrex::Real z)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    Tensor::Rank1 cartesian_v_L{0., 0., 0.};

    // derivatives for jacobian matrix - drdx etc
    Tensor::Rank2 jac = spherical_jacobian(x, y, z);

    // transform the vector to cartesian coords
    FOR (i)
    {
        FOR (j)
        {
            cartesian_v_L(i) += spherical_v_L(j) * jac(j, i);
        }
    }
    return cartesian_v_L;
}

// Convert a vector (with one upper index) in cartesian coords to spherical
// coords

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
Tensor::Rank1 cartesian_to_spherical_U(const Tensor::Rank1 &cartesian_v_U,
                                       amrex::Real x, amrex::Real y,
                                       amrex::Real z)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    Tensor::Rank1 spherical_v_U{0., 0., 0.};

    // derivatives for jacobian matrix - drdx etc
    Tensor::Rank2 jac = spherical_jacobian(x, y, z);

    // transform the vector to cartesian coords
    FOR (i)
    {
        FOR (j)
        {
            spherical_v_U(i) += jac(i, j) * cartesian_v_U(j);
        }
    }
    return spherical_v_U;
}

// Convert a vector (with one lower index) in cartesian coords to spherical
// coords

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
Tensor::Rank1 cartesian_to_spherical_L(const Tensor::Rank1 &cartesian_v_L,
                                       amrex::Real x, amrex::Real y,
                                       amrex::Real z)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    Tensor::Rank1 spherical_v_L{0., 0., 0.};

    // derivatives for inverse jacobian matrix - drdx etc
    Tensor::Rank2 inv_jac = inverse_spherical_jacobian(x, y, z);

    // transform the vector to cartesian coords
    FOR (i)
    {
        FOR (j)
        {
            spherical_v_L(i) += cartesian_v_L(j) * inv_jac(j, i);
        }
    }
    return spherical_v_L;
}

// The area element of a sphere
amrex::Real area_element_sphere(const Tensor::Sym12Rank2 &spherical_g)
{
    return std::sqrt(spherical_g(1, 1) * spherical_g(2, 2) -
                     spherical_g(1, 2) * spherical_g(2, 1));
}

} // namespace CoordinateTransformations
#endif /* COORDINATETRANSFORMATIONS_HPP_ */
