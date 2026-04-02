/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(TWOPUNCTURESINITIALDATA_HPP_)
#error "This file should only be included through TwoPuncturesInitialData.hpp"
#endif

#ifndef TWOPUNCTURESINITIALDATA_IMPL_HPP_
#define TWOPUNCTURESINITIALDATA_IMPL_HPP_

void TwoPuncturesInitialData::compute(Cell<double> current_cell) const
{
    Vars<amrex::Real> vars;
    // Set only the non-zero components explicitly below
    VarsTools::assign(vars, 0.);

    Coordinates<amrex::Real> coords(current_cell, m_dx, m_center);
    TensorArray::Rank2 h_phys, K_tensor;
    TensorArray::Rank1 shift, Z3;
    double lapse, Theta;

    interpolate_tp_vars(coords, h_phys, K_tensor, lapse, shift, Theta, Z3);

    using namespace TensorAlgebra;
    // analytically set Bowen-York properties below (e.g. conformal flatness,
    // tracefree K)

    // metric variables
    vars.chi = pow(compute_determinant_sym(h_phys), -1. / 3.);
    FOR (i)
    {
        // Bowen-York data is conformally flat
        vars.h(i, i) = 1.0;
    }

    // extrinsic curvature
    FOR (i, j)
    {
        vars.A(i, j) = vars.chi * K_tensor(i, j);
    }
    // conformal flatness means h_UU = h
    make_trace_free(vars.A, vars.h, vars.h);

    // gauge
    vars.lapse = lapse;

    current_cell.store_vars(vars);
}

void TwoPuncturesInitialData::interpolate_tp_vars(
    const Coordinates<amrex::Real> &coords, TensorArray::Rank2 &out_h_phys,
    TensorArray::Rank2 &out_K_tensor, amrex::Real &out_lapse,
    TensorArray::Rank1 &out_shift, amrex::Real &out_Theta,
    TensorArray::Rank1 &out_Z3) const
{
    amrex::Real coords_array[AMREX_SPACEDIM];
    coords_array[0] = coords.x;
    coords_array[1] = coords.y;
    coords_array[2] = coords.z;

    using namespace TP::Z4VectorShortcuts;
    amrex::Real TP_state[Qlen];
    m_two_punctures.Interpolate(coords_array, TP_state);

    // metric
    out_h_phys(0, 0) = TP_state[g11];
    out_h_phys(0, 1) = out_h_phys(1, 0) = TP_state[g12];
    out_h_phys(0, 2) = out_h_phys(2, 0) = TP_state[g13];
    out_h_phys(1, 1)                    = TP_state[g22];
    out_h_phys(1, 2) = out_h_phys(2, 1) = TP_state[g23];
    out_h_phys(2, 2)                    = TP_state[g33];

    // extrinsic curvature
    out_K_tensor(0, 0) = TP_state[K11];
    out_K_tensor(0, 1) = out_K_tensor(1, 0) = TP_state[K12];
    out_K_tensor(0, 2) = out_K_tensor(2, 0) = TP_state[K13];
    out_K_tensor(1, 1)                      = TP_state[K22];
    out_K_tensor(1, 2) = out_K_tensor(2, 1) = TP_state[K23];
    out_K_tensor(2, 2)                      = TP_state[K33];

    // Z4 vector
    out_Z3(0) = TP_state[Z1];
    out_Z3(1) = TP_state[Z2];
    out_Z3(2) = TP_state[Z3];
    out_Theta = TP_state[Theta];

    // gauge
    out_lapse    = TP_state[lapse];
    out_shift(0) = TP_state[shift1];
    out_shift(1) = TP_state[shift2];
    out_shift(2) = TP_state[shift3];
}

#endif /* TWOPUNCTURESINITIALDATA_IMPL_HPP_ */
