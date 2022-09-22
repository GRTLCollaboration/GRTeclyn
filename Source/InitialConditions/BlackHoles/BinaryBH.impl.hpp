/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(BINARYBH_HPP_)
#error "This file should only be included through BinaryBH.hpp"
#endif

#ifndef BINARYBH_IMPL_HPP_
#define BINARYBH_IMPL_HPP_

#include "BSSNVars.hpp"
#include "BinaryBH.hpp"
#include "VarsTools.hpp"
#include "simd.hpp"

template <class data_t> void BinaryBH::compute(Cell<data_t> current_cell) const
{
    BSSNVars::VarsWithGauge<data_t> vars;
    VarsTools::assign(vars,
                      0.); // Set only the non-zero components explicitly below
    Coordinates<data_t> coords(current_cell, m_dx);

    vars.chi = compute_chi(coords);

    // Conformal metric is flat
    FOR(i) vars.h[i][i] = 1.;

    vars.A = compute_A(vars.chi, coords);

    switch (m_initial_lapse)
    {
    case Lapse::ONE:
        vars.lapse = 1.;
        break;
    case Lapse::PRE_COLLAPSED:
        vars.lapse = sqrt(vars.chi);
        break;
    case Lapse::CHI:
        vars.lapse = vars.chi;
        break;
    default:
        amrex::Abort("BinaryBH::Supplied initial lapse not supported.");
    }

    current_cell.store_vars(vars);
}

template <class data_t>
data_t BinaryBH::compute_chi(Coordinates<data_t> coords) const
{
    const data_t psi =
        1. + bh1.psi_minus_one(coords) + bh2.psi_minus_one(coords);
    return pow(psi, -4);
}

template <class data_t>
Tensor<2, data_t> BinaryBH::compute_A(data_t chi,
                                      Coordinates<data_t> coords) const
{

    Tensor<2, data_t> Aij1 = bh1.Aij(coords);
    Tensor<2, data_t> Aij2 = bh2.Aij(coords);
    Tensor<2, data_t> out;

    // Aij(CCZ4) = psi^(-6) * Aij(Baumgarte&Shapiro book)
    FOR(i, j) out[i][j] = pow(chi, 3 / 2.) * (Aij1[i][j] + Aij2[i][j]);

    return out;
}

template <class data_t>
AMREX_GPU_DEVICE
void BinaryBH::init_data (int i, int j, int k,
                          amrex::CellData<data_t> const& cell) const
{
    BSSNVars::VarsWithGauge<data_t> vars;
    VarsTools::assign(vars,
                      0.); // Set only the non-zero components explicitly below
    Coordinates<data_t> coords(amrex::IntVect(i,j,k), m_dx);

    vars.chi = compute_chi(coords);

    // Conformal metric is flat
    FOR(i) vars.h[i][i] = 1.;

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
        amrex::Abort("BinaryBH::Supplied initial lapse not supported.");
    }

    store_vars(cell, vars);
}

#endif /* BINARYBH_IMPL_HPP_ */
