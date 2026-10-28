/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef COORDINATETRANSFORMATIONS_HPP_
#define COORDINATETRANSFORMATIONS_HPP_

#include "DimensionDefinitions.hpp"
#include "TensorAlgebra.hpp"

#include <cmath>

#include <AMReX_Math.H>

using namespace amrex::literals;

namespace CoordinateTransformations
{

// Jacobian transformation matrix

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor::Rank2
spherical_jacobian(const amrex::Real x, const amrex::Real y,
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
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor::Rank2
inverse_spherical_jacobian(const amrex::Real x, const amrex::Real y,
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

// transform a vector according to a given jacobian
// NOLINTBEGIN(bugprone-easily-swappable-parameters)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor::Rank1
transform_vector(const Tensor::Rank1 &vec_U, const Tensor::Rank2 &jacobian)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    Tensor::Rank1 transformed_U{0.0, 0.0, 0.0};
    FOR (i)
    {
        FOR (j)
        {
            transformed_U(i) += jacobian(i, j) * vec_U(j);
        }
    }
    return transformed_U;
}

// transform a covector according to a given jacobian
// NOLINTBEGIN(bugprone-easily-swappable-parameters)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor::Rank1
transform_covector(const Tensor::Rank1 &vec_L, const Tensor::Rank2 &jacobian)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    Tensor::Rank1 transformed_L{0.0, 0.0, 0.0};
    FOR (i)
    {
        FOR (j)
        {
            transformed_L(i) += vec_L(j) * jacobian(j, i);
        }
    }
    return transformed_L;
}

// transform a full tensor with two raised indices
// NOLINTBEGIN(bugprone-easily-swappable-parameters)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor::Rank2
transform_tensor_UU(const Tensor::Rank2 &tensor_UU,
                    const Tensor::Rank2 &jacobian)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    Tensor::Rank2 transformed_UU{0.0};
    FOR (i, j)
    {
        FOR (k, l)
        {
            transformed_UU(i, j) +=
                jacobian(i, k) * jacobian(j, l) * tensor_UU(k, l);
        }
    }
    return transformed_UU;
}

// transform a full tensor with two lowered indices
// NOLINTBEGIN(bugprone-easily-swappable-parameters)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor::Rank2
transform_tensor_LL(const Tensor::Rank2 &tensor_LL,
                    const Tensor::Rank2 &jacobian)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    Tensor::Rank2 transformed_LL{0.0};
    FOR (i, j)
    {
        FOR (k, l)
        {
            transformed_LL(i, j) +=
                tensor_LL(k, l) * jacobian(k, i) * jacobian(l, j);
        }
    }
    return transformed_LL;
}

// transform a *symmetric* tensor with two raised indices
// NOLINTBEGIN(bugprone-easily-swappable-parameters)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor::Sym12Rank2
transform_sym_tensor_UU(const Tensor::Sym12Rank2 &tensor_UU,
                        const Tensor::Rank2 &jacobian)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    Tensor::Sym12Rank2 transformed_UU{0.0};
    FOR2_SYM(i, j)
    {
        FOR (k, l)
        {
            transformed_UU(i, j) +=
                jacobian(i, k) * jacobian(j, l) * tensor_UU(k, l);
        }
    }
    return transformed_UU;
}

// transform a *symmetric* tensor with two lowered indices
// NOLINTBEGIN(bugprone-easily-swappable-parameters)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor::Sym12Rank2
transform_sym_tensor_LL(const Tensor::Sym12Rank2 &tensor_LL,
                        const Tensor::Rank2 &jacobian)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    Tensor::Sym12Rank2 transformed_LL{0.0};
    FOR2_SYM(i, j)
    {
        FOR (k, l)
        {
            transformed_LL(i, j) +=
                tensor_LL(k, l) * jacobian(k, i) * jacobian(l, j);
        }
    }
    return transformed_LL;
}

// Convert a Tensor (with two lower indices) in spherical coords to cartesian
// coords

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor::Sym12Rank2
spherical_to_cartesian_LL(const Tensor::Sym12Rank2 &spherical_g,
                          const amrex::Real x, const amrex::Real y,
                          const amrex::Real z)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    return transform_sym_tensor_LL(spherical_g, spherical_jacobian(x, y, z));
}

// Convert a Tensor (with two upper indices) in spherical coords to cartesian
// coords

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor::Sym12Rank2
spherical_to_cartesian_UU(const Tensor::Sym12Rank2 &spherical_g_UU,
                          const amrex::Real x, const amrex::Real y,
                          const amrex::Real z)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    return transform_sym_tensor_UU(spherical_g_UU,
                                   inverse_spherical_jacobian(x, y, z));
}

// Convert a Tensor (with two lower indices) in cartesian coords to spherical
// coords

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor::Sym12Rank2
cartesian_to_spherical_LL(const Tensor::Sym12Rank2 &cartesian_g,
                          const amrex::Real x, const amrex::Real y,
                          const amrex::Real z)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    return transform_sym_tensor_LL(cartesian_g,
                                   inverse_spherical_jacobian(x, y, z));
}

// Convert a Tensor (with two upper indices) in cartesian coords to spherical
// coords

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor::Sym12Rank2
cartesian_to_spherical_UU(const Tensor::Sym12Rank2 &cartesian_g_UU,
                          amrex::Real x, amrex::Real y, amrex::Real z)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    return transform_sym_tensor_UU(cartesian_g_UU, spherical_jacobian(x, y, z));
}

// Convert a vector (with one upper index) in spherical coords to cartesian
// coords

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor::Rank1
spherical_to_cartesian_U(const Tensor::Rank1 &spherical_v_U, amrex::Real x,
                         amrex::Real y, amrex::Real z)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    return transform_vector(spherical_v_U, inverse_spherical_jacobian(x, y, z));
}

// Convert a vector (with one lower index) in spherical coords to cartesian
// coords

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor::Rank1
spherical_to_cartesian_L(const Tensor::Rank1 &spherical_v_L, amrex::Real x,
                         amrex::Real y, amrex::Real z)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    return transform_covector(spherical_v_L, spherical_jacobian(x, y, z));
}

// Convert a vector (with one upper index) in cartesian coords to spherical
// coords

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor::Rank1
cartesian_to_spherical_U(const Tensor::Rank1 &cartesian_v_U, amrex::Real x,
                         amrex::Real y, amrex::Real z)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    return transform_vector(cartesian_v_U, spherical_jacobian(x, y, z));
}

// Convert a vector (with one lower index) in cartesian coords to spherical
// coords

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor::Rank1
cartesian_to_spherical_L(const Tensor::Rank1 &cartesian_v_L, amrex::Real x,
                         amrex::Real y, amrex::Real z)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    return transform_covector(cartesian_v_L, inverse_spherical_jacobian(x, y, z));
}

// The area element of a sphere
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
area_element_sphere(const Tensor::Sym12Rank2 &spherical_g)
{
    return std::sqrt(spherical_g(1, 1) * spherical_g(2, 2) -
                     spherical_g(1, 2) * spherical_g(2, 1));
}

// cartesian coordinates rotation matrix
// axis must be normalized to 1
// NOLINTBEGIN(bugprone-easily-swappable-parameters)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor::Rank2
rotation_matrix(const Tensor::Rank1 &axis, amrex::Real cos_angle)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    amrex::Real sine =
        std::sqrt(std::max(0._rt, 1. - cos_angle * cos_angle));
    amrex::Real one_minus_cos = 1. - cos_angle;

    Tensor::Rank2 rotation;
    rotation(0, 0) = axis(0) * axis(0) * one_minus_cos + cos_angle;
    rotation(0, 1) = axis(1) * axis(0) * one_minus_cos - axis(2) * sine;
    rotation(0, 2) = axis(2) * axis(0) * one_minus_cos + axis(1) * sine;
    rotation(1, 0) = axis(0) * axis(1) * one_minus_cos + axis(2) * sine;
    rotation(1, 1) = axis(1) * axis(1) * one_minus_cos + cos_angle;
    rotation(1, 2) = axis(2) * axis(1) * one_minus_cos - axis(0) * sine;
    rotation(2, 0) = axis(0) * axis(2) * one_minus_cos - axis(1) * sine;
    rotation(2, 1) = axis(1) * axis(2) * one_minus_cos + axis(0) * sine;
    rotation(2, 2) = axis(2) * axis(2) * one_minus_cos + cos_angle;

    return rotation;
}

// cartesian coordinates rotation matrix from origin vector to destination
// vector
// NOLINTBEGIN(bugprone-easily-swappable-parameters)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor::Rank2
rotation_matrix(const Tensor::Rank1 &origin, const Tensor::Rank1 &destination)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    // Calculate axis of rotation:
    Tensor::Rank1 axis = {
        origin(1) * destination(2) - origin(2) * destination(1),
        origin(2) * destination(0) - origin(0) * destination(2),
        origin(0) * destination(1) - origin(1) * destination(0)};

    // Calculate angle of rotation
    amrex::Real norm_origin = std::sqrt(
        origin(0) * origin(0) + origin(1) * origin(1) + origin(2) * origin(2));
    amrex::Real norm_dest = std::sqrt(destination(0) * destination(0) +
                                      destination(1) * destination(1) +
                                      destination(2) * destination(2));

    amrex::Real cos_angle  = origin(0) * destination(0) +
                             origin(1) * destination(1) +
                             origin(2) * destination(2);
    cos_angle             /= (norm_origin * norm_dest);

    // origin and destination are (anti-)parallel so the axis is undefined
    constexpr amrex::Real eps = 1.e-13;
    amrex::Real axis_norm2 =
        axis(0) * axis(0) + axis(1) * axis(1) + axis(2) * axis(2);
    if (axis_norm2 <
        eps * eps * norm_origin * norm_origin * norm_dest * norm_dest)
    {
        if (cos_angle > 0.)
        {
            return rotation_matrix({0., 0., 0.}, 1.);
        }
        // anti-parallel: rotate by pi about any axis orthogonal to origin,
        // e.g. origin x e_k with e_k the cartesian axis least aligned with it
        int k = 0;
        for (int i = 1; i < AMREX_SPACEDIM; ++i)
        {
            if (std::abs(origin(i)) < std::abs(origin(k)))
            {
                k = i;
            }
        }
        Tensor::Rank1 e_k = {0., 0., 0.};
        e_k(k)            = 1.;
        axis              = {origin(1) * e_k(2) - origin(2) * e_k(1),
                             origin(2) * e_k(0) - origin(0) * e_k(2),
                             origin(0) * e_k(1) - origin(1) * e_k(0)};
        axis_norm2 = axis(0) * axis(0) + axis(1) * axis(1) + axis(2) * axis(2);
        cos_angle  = -1.;
    }

    amrex::Real norm_inv = amrex::Math::rsqrt(axis_norm2);

    axis(0) *= norm_inv;
    axis(1) *= norm_inv;
    axis(2) *= norm_inv;

    return rotation_matrix(axis, cos_angle);
}

} // namespace CoordinateTransformations
#endif /* COORDINATETRANSFORMATIONS_HPP_ */
