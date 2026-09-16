/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(KERRBHINITIALDATA_HPP_)
#error "This file should only be included through KerrBHInitialData.hpp"
#endif

#ifndef KERRBHINITIALDATA_IMPL_HPP_
#define KERRBHINITIALDATA_IMPL_HPP_

#include "CCZ4Vars.hpp"
#include "DimensionDefinitions.hpp"
#include "GRParmParse.hpp"
#include <cmath>

inline void KerrBHInitialData::params_t::check_params()
{
    GRParmParse kerr_pp("kerr");
    amrex::Real check_mass{};
    kerr_pp.get("mass", check_mass);
    
    amrex::Real check_spin{};
    kerr_pp.get("spin", check_spin);

    // Validate physics before the run starts
    if (std::abs(check_spin) > check_mass)
    {
        amrex::Abort("The spin parameter must satisfy |a| <= M");
    }

    // Ensure we don't crash if geometry.center is missing from params.txt
    GRParmParse geom_pp("geometry");
    std::array<amrex::Real, AMREX_SPACEDIM> check_center{};
    geom_pp.queryAdd("center", check_center);
}

inline void KerrBHInitialData::params_t::fill_params()
{
    GRParmParse kerr_pp("kerr");
    kerr_pp.get("mass", mass);
    kerr_pp.get("spin", spin);
    kerr_pp.queryAdd("spin_direction", spin_direction);

    GRParmParse geom_pp("geometry");
    geom_pp.queryAdd("center", center);

    std::array<amrex::Real, AMREX_SPACEDIM> offset{};
    kerr_pp.queryAdd("offset", offset);
    FOR (idir)
    {
        center[idir] += offset[idir];
    }
}

AMREX_FORCE_INLINE
KerrBHInitialData::KerrBHInitialData(amrex::Real a_dx) : m_dx(a_dx)
{
    m_params.fill_params();
    
    if (std::abs(m_params.spin) > m_params.mass)
    {
        amrex::Abort("The spin parameter must satisfy |a| <= M");
    }

    // define the rotation matrix needed to transform Cartesian
    // coordinates into the coordinates of the spin direction
    Tensor::Rank1 z_dir = {0., 0., 1.};
    Tensor::Rank1 spin_dir = {m_params.spin_direction[0],
                              m_params.spin_direction[1],
                              m_params.spin_direction[2]};
    m_R = CoordinateTransformations::rotation_matrix(spin_dir, z_dir);
}

// Unpack symmetric 6-element array into a full 3x3 matrix for TensorAlgebra
AMREX_FORCE_INLINE AMREX_GPU_DEVICE Tensor::Rank2
KerrBHInitialData::to_rank2(const Tensor::Sym12Rank2& sym) const
{
    Tensor::Rank2 r2;
    FOR2_SYM(i, j) 
    { 
        r2(i, j) = r2(j, i) = sym(i, j);
    }
    return r2;
}

// Computes semi-isotropic Kerr solution as detailed in Liu, Etienne and Shapiro
// 2010, arxiv gr-qc/1001.4077
AMREX_FORCE_INLINE AMREX_GPU_DEVICE void
KerrBHInitialData::operator()(int ix, int iy, int iz, 
                              const amrex::Array4<amrex::Real> &state) const
{
    using namespace CoordinateTransformations;
    using namespace TensorAlgebra;

    // set up vars for the metric and extrinsic curvature, shift and lapse in
    // spherical coords

    Tensor::Sym12Rank2 spherical_g{0.};
    Tensor::Sym12Rank2 spherical_K{0.};
    Tensor::Rank1 spherical_shift{0., 0., 0.};
    amrex::Real kerr_lapse;

    // The cartesian variables and coords
    const amrex::CellData<amrex::Real> &state_cell_data = state.cellData(ix, iy, iz);
    Coordinates coords(amrex::IntVect(ix, iy, iz), m_dx);

    Tensor::Rank1 xyz = {coords.x - m_params.center[0],
                         coords.y - m_params.center[1],
                         coords.z - m_params.center[2]};

    xyz = transform_vector(xyz, m_R);

    // Compute the components in spherical coords as per 1401.1548
    compute_kerr(spherical_g, spherical_K, spherical_shift, kerr_lapse, xyz);

    // work out where we are on the grid
    amrex::Real x = xyz(0);
    amrex::Real y = xyz(1);
    amrex::Real z = xyz(2);

    // Convert spherical components to cartesian components using coordinate
    // transforms
    auto cartesian_h_sym = spherical_to_cartesian_LL(spherical_g, x, y, z);
    auto cartesian_K_sym = spherical_to_cartesian_LL(spherical_K, x, y, z);
    auto cartesian_shift = spherical_to_cartesian_U(spherical_shift, x, y, z);

    // Rotate symmetric tensors back to original coordinates
    auto h_sym = transform_sym_tensor_LL(cartesian_h_sym, m_R);
    auto A_sym = transform_sym_tensor_LL(cartesian_K_sym, m_R);
    auto shift = transform_vector(cartesian_shift, m_R);

    // Unpack Sym12Rank2 into Rank2 for trace and determinant operations
    Tensor::Rank2 h = to_rank2(h_sym);
    Tensor::Rank2 A = to_rank2(A_sym);

    // Convert to BSSN vars
    amrex::Real deth = compute_determinant(h);
    auto h_UU_sym    = compute_inverse_sym(h_sym);
    
    // Unpack Inverse to Rank2
    Tensor::Rank2 h_UU = to_rank2(h_UU_sym);

    amrex::Real chi  = pow(deth, -1. / 3.);

    // transform extrinsic curvature into A and TrK - note h is still non
    // conformal version which is what we need here
    amrex::Real K = compute_trace(A, h_UU);
    make_trace_free(A, h, h_UU);

    // Make conformal
    FOR2_SYM(i, j)
    {
        h(i, j) *= chi;
        A(i, j) *= chi;
    }

    // use a pre collapsed lapse, could also use analytic one
    // amrex::Real lapse = kerr_lapse;
    amrex::Real lapse = pow(chi, 0.5);

    // Populate the variables on the grid
    // NB We still need to set Gamma^i which is NON ZERO
    // but we do this via a separate class/compute function
    // as we need the gradients of the metric which are not yet available
    state_cell_data[c_chi]   = chi;
    state_cell_data[c_K]     = K;
    state_cell_data[c_lapse] = lapse;

    FOR2_SYM(i, j)
    {
        state_cell_data[sym_var_idx(c_h11, i, j)] = h(i, j);
        state_cell_data[sym_var_idx(c_A11, i, j)] = A(i, j);
    }

    FOR(i)
    {
        state_cell_data[c_shift1 + i] = shift(i);
    }
}

AMREX_FORCE_INLINE AMREX_GPU_DEVICE void 
KerrBHInitialData::compute_kerr(Tensor::Sym12Rank2 &spherical_g,
                                Tensor::Sym12Rank2 &spherical_K,
                                Tensor::Rank1 &spherical_shift,
                                amrex::Real &kerr_lapse,
                                const Tensor::Rank1 &coords) const
{
    // Kerr black hole params - mass M and spin a
    amrex::Real M = m_params.mass;
    amrex::Real a = m_params.spin;

    // work out where we are on the grid
    amrex::Real x = coords(0);
    amrex::Real y = coords(1);
    amrex::Real z = coords(2);

    // the radius, subject to a floor
    amrex::Real r2_raw = x * x + y * y + z * z;
    amrex::Real r = std::max(sqrt(r2_raw), 1e-6);
    amrex::Real r2 = r * r;

    // the radius in xy plane, subject to a floor
    amrex::Real rho2 = std::max(x * x + y * y, 1e-12);
    amrex::Real rho  = sqrt(rho2);

    // calculate useful position quantities
    amrex::Real cos_theta  = z / r;
    amrex::Real sin_theta  = rho / r;
    amrex::Real cos_theta2 = cos_theta * cos_theta;
    amrex::Real sin_theta2 = sin_theta * sin_theta;

    // calculate useful metric quantities
    amrex::Real r_plus  = M + sqrt(M * M - a * a);
    amrex::Real r_minus = M - sqrt(M * M - a * a);

    // The Boyer-Lindquist coordinate
    amrex::Real r_BL = r * pow(1.0 + 0.25 * r_plus / r, 2.0);

    // Other useful quantities per 1001.4077
    amrex::Real Sigma = r_BL * r_BL + a * a * cos_theta2;
    amrex::Real Delta = r_BL * r_BL - 2.0 * M * r_BL + a * a;
    // In the paper this is just 'A', but not to be confused with A_ij
    amrex::Real AA = pow(r_BL * r_BL + a * a, 2.0) - Delta * a * a * sin_theta2;
    // The rr component of the conformal spatial matric
    amrex::Real gamma_rr =
        Sigma * pow(r + 0.25 * r_plus, 2.0) / (r * r2 * (r_BL - r_minus));

    // Zero initialize the symmetric rank 2 arrays safely
    FOR2_SYM (i, j) 
    { 
        spherical_g(i, j) = 0.0; 
        spherical_K(i, j) = 0.0;
    }

    // Metric in semi isotropic Kerr-Schild coordinates, r, theta (t or th), phi
    // (p)
    spherical_g(0, 0) = gamma_rr;                // gamma_rr
    spherical_g(1, 1) = Sigma;                   // gamma_tt
    spherical_g(2, 2) = AA / Sigma * sin_theta2; // gamma_pp

    // Extrinsic curvature
    // set non zero elements of Krtp - K_rp, K_tp
    spherical_K(0, 2) =
        a * M * sin_theta2 / (Sigma * sqrt(AA * Sigma)) *
        (3.0 * pow(r_BL, 4.0) + 2 * a * a * r_BL * r_BL - pow(a, 4.0) -
         a * a * (r_BL * r_BL - a * a) * sin_theta2) *
        (1.0 + 0.25 * r_plus / r) / sqrt(r * r_BL - r * r_minus);
    
    spherical_K(1, 2) = 
        -2.0 * pow(a, 3.0) * M * r_BL * cos_theta * sin_theta *
        sin_theta2 / (Sigma * sqrt(AA * Sigma)) *
        (r - 0.25 * r_plus) * sqrt(r_BL / r - r_minus / r);

    // set the analytic lapse
    kerr_lapse = sqrt(Delta * Sigma / AA);

    // set the shift (only the phi component is non zero)
    spherical_shift(0) = 0.0;
    spherical_shift(1) = 0.0;
    spherical_shift(2) = -2.0 * M * a * r_BL / AA;
}

#endif /* KERRBHINITIALDATA_IMPL_HPP_ */
