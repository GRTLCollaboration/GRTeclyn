/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(BINARYBHINITIALDATA_HPP_)
#error "This file should only be included through BinaryBHInitialData.hpp"
#endif

#ifndef BINARYBHINITIALDATA_IMPL_HPP_
#define BINARYBHINITIALDATA_IMPL_HPP_

#include "BSSNVars.hpp"
#include "VarsTools.hpp"

// Constructor
// NOLINTBEGIN(bugprone-easily-swappable-parameters)
AMREX_FORCE_INLINE BinaryBHInitialData::BinaryBHInitialData(
    BoostedBHInitialData::params_t a_bh1_params,
    BoostedBHInitialData::params_t a_bh2_params, double a_dx,
    int a_initial_lapse)
    : m_dx(a_dx), bh1(a_bh1_params), bh2(a_bh2_params),
      m_initial_lapse(a_initial_lapse)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
}

[[nodiscard]] AMREX_FORCE_INLINE AMREX_GPU_DEVICE amrex::Real
BinaryBHInitialData::compute_chi(Coordinates coords) const
{
    const amrex::Real psi =
        1. + bh1.psi_minus_one(coords) + bh2.psi_minus_one(coords);
    return pow(psi, -4);
}

[[nodiscard]] AMREX_FORCE_INLINE AMREX_GPU_DEVICE Tensor<2, amrex::Real>
BinaryBHInitialData::compute_A(amrex::Real chi, Coordinates coords) const
{

    Tensor<2, amrex::Real> Aij1 = bh1.Aij(coords);
    Tensor<2, amrex::Real> Aij2 = bh2.Aij(coords);
    Tensor<2, amrex::Real> out;

    // Aij(CCZ4) = psi^(-6) * Aij(Baumgarte&Shapiro book)
    FOR (i, j)
        out[i][j] = pow(chi, 3 / 2.) * (Aij1[i][j] + Aij2[i][j]);

    return out;
}

AMREX_FORCE_INLINE
AMREX_GPU_DEVICE // or AMREX_GPU_HOST_DEVICE depending on what's needed
    void
    BinaryBHInitialData::init_data(
        int i, int j, int k, const amrex::CellData<amrex::Real> &cell) const
{
    // TODO: Remove this once BSSNVars de-data_t-ed
    BSSNVars::VarsWithGauge<amrex::Real> vars;
    VarsTools::assign(vars,
                      0.); // Set only the non-zero components explicitly below
    Coordinates coords(amrex::IntVect(i, j, k), m_dx);

    vars.chi = compute_chi(coords);

    // Conformal metric is flat
    FOR (ii)
        vars.h[ii][ii] = 1.;

    vars.A = compute_A(vars.chi, coords);

    switch (m_initial_lapse)
    {
    case Lapse::ONE:
        vars.lapse = 1.;
        break;
    case Lapse::PRE_COLLAPSED:
        vars.lapse = std::sqrt(vars.chi);
        break;
    case Lapse::CHI:
        vars.lapse = vars.chi;
        break;
    default:
        amrex::Abort(
            "BinaryBHInitialData::Supplied initial lapse not supported.");
    }

    store_vars(cell, vars);
}

#endif /* BINARYBHINITIALDATA_IMPL_HPP_ */
