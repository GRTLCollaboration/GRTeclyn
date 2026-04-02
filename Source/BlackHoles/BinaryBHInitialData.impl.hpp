/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(BINARYBHINITIALDATA_HPP_)
#error "This file should only be included through BinaryBHInitialData.hpp"
#endif

#ifndef BINARYBHINITIALDATA_IMPL_HPP_
#define BINARYBHINITIALDATA_IMPL_HPP_

#include "BinaryBHInitialData.hpp"
#include "BoostedBHInitialData.hpp"
#include "CCZ4Vars.hpp"
#include "TensorAlgebra.hpp"

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

[[nodiscard]] AMREX_FORCE_INLINE AMREX_GPU_DEVICE TensorArray::Rank2
BinaryBHInitialData::compute_A(amrex::Real chi, Coordinates coords) const
{

    TensorArray::Rank2 Aij1 = bh1.Aij(coords);
    TensorArray::Rank2 Aij2 = bh2.Aij(coords);
    TensorArray::Rank2 Aij_total;

    // Aij(CCZ4) = psi^(-6) * Aij(Baumgarte&Shapiro book)
    FOR (i, j)
        Aij_total(i, j) = pow(chi, 3 / 2.) * (Aij1(i, j) + Aij2(i, j));

    return Aij_total;
}

AMREX_FORCE_INLINE
AMREX_GPU_DEVICE // or AMREX_GPU_HOST_DEVICE depending on what's needed
    void
    BinaryBHInitialData::operator()(
        int ix, int iy, int iz, const amrex::Array4<amrex::Real> &state) const
{

    const amrex::CellData<amrex::Real> &state_cell_data =
        state.cellData(ix, iy, iz);
    Coordinates coords(amrex::IntVect(ix, iy, iz), m_dx);

    // Assign non zero values
    amrex::Real chi        = compute_chi(coords);
    state_cell_data[c_chi] = chi;

    FOR2_SYM(i, j)
    {
        state_cell_data[VAR_IDX(c_h11, i, j)] = TensorAlgebra::delta(i, j);
    }

    auto total_A_LL = compute_A(chi, coords);
    FOR2_SYM(i, j) { state_cell_data[VAR_IDX(c_A11, i, j)] = total_A_LL(i, j); }

    switch (m_initial_lapse)
    {
    case Lapse::ONE:
        state_cell_data[c_lapse] = 1.0;
        break;
    case Lapse::PRE_COLLAPSED:
        state_cell_data[c_lapse] = std::sqrt(chi);
        break;
    case Lapse::CHI:
        state_cell_data[c_lapse] = chi;
        break;
    default:
        amrex::Abort(
            "BinaryBHInitialData::Supplied initial lapse not supported.");
    }
}

#endif /* BINARYBHINITIALDATA_IMPL_HPP_ */
