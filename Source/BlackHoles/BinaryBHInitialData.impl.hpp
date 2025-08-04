/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#include "BinaryBHInitialData.hpp"
#include "BoostedBHInitialData.hpp"
#include "CCZ4Vars2.hpp"
#include "TensorAlgebra2.hpp"
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
    Tensor<2, amrex::Real> Aij_total;

    // Aij(CCZ4) = psi^(-6) * Aij(Baumgarte&Shapiro book)
    FOR (i, j)
        Aij_total[i][j] = pow(chi, 3 / 2.) * (Aij1[i][j] + Aij2[i][j]);

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
    CCZ4Vars2 vars(state_cell_data);
    Coordinates coords(amrex::IntVect(ix, iy, iz), m_dx);

    // Assign non zero values
    amrex::Real chi = compute_chi(coords);
    vars.store_chi(chi);

    Tensor<2, amrex::Real> h_LL;
    FOR2 (i, j)
    {
        h_LL[i][j] = TensorAlgebra2::delta(i, j);
    }
    vars.store_h(h_LL);

    Tensor<2, amrex::Real> total_A_LL = compute_A(chi, coords);
    vars.store_A(total_A_LL);

    switch (m_initial_lapse)
    {
    case Lapse::ONE:
        vars.store_lapse(1.0);
        break;
    case Lapse::PRE_COLLAPSED:
        vars.store_lapse(std::sqrt(chi));
        break;
    case Lapse::CHI:
        vars.store_lapse(chi);
        break;
    default:
        amrex::Abort(
            "BinaryBHInitialData::Supplied initial lapse not supported.");
    }
}

#endif /* BINARYBHINITIALDATA_IMPL_HPP_ */
