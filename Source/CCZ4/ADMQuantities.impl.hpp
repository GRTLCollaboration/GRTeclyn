/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(ADMQUANTITIES_HPP_)
#error "This file should only be included through ADMQuantities.hpp"
#endif

#ifndef ADMQUANTITIES_IMPL_HPP_
#define ADMQUANTITIES_IMPL_HPP_

#include "ADMQuantities.hpp"

// AMReX includes
#include <AMReX_AmrLevel.H>

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
AMREX_FORCE_INLINE
ADMQuantities::ADMQuantities(amrex::Real dx, int a_c_Madm) 
: m_deriv(dx), m_c_Madm(a_c_Madm),m_center{}
// m_G_Newton(a_G_Newton),
// m_c_Madm(a_c_Madm), m_c_Jadm(a_c_Jadm), m_dir(Z)
{
    m_G_Newton = 1.0;
    m_dx =dx;
    // if (a_G_Newton == 0) || a_G_Newton not supplied
    // {
    //     m_G_Newton = a_G_Newton;
    // }
}
// NOLINTEND(bugprone-easily-swappable-parameters)

AMREX_GPU_DEVICE void
ADMQuantities::operator()(int ix, int iy, int iz,
                        const amrex::Array4<amrex::Real> &adm_state,
                        const amrex::Array4<const amrex::Real> &state) const
{
    // We do not want to amend the cell data values, so use const CCZ4Vars
    // const amrex::CellData<const amrex::Real> &state_cell_data =state.cellData(ix, iy, iz);
    amrex::IntVect cell(AMREX_D_DECL(ix, iy, iz));

    const amrex::CellData<const amrex::Real> &state_cell_data =
        state.cellData(ix, iy, iz);
    CCZ4Vars vars(state_cell_data);

    const auto h_UU  = CCZ4Geometry::compute_inverse_metric(vars);

    // we need d1 chi, K, h, A... this just gets all of them
    Tensor::Rank1 d1_chi    = m_deriv.d1_scalar(ix, iy, iz, state, c_chi);
    // Tensor::Rank2 d1_Gamma  = m_deriv.d1_vector(ix, iy, iz, state, c_Gamma1);
    // Tensor::Rank1 d1_K      = m_deriv.d1_scalar(ix, iy, iz, state, c_K);
    // Tensor::Sym12Rank3 d1_A = m_deriv.d1_sym_tensor(ix, iy, iz, state, c_A11);
    Tensor::Sym12Rank3 d1_h = m_deriv.d1_sym_tensor(ix, iy, iz, state, c_h11);

    // Surface element for integration
    Coordinates coords(cell, m_dx, m_center);
    Tensor::Rank1 x = {coords.x, coords.y, coords.z};
    Tensor::Rank1 dS_U = x;

    amrex::Real dS_norm{};
    FOR(i, j) { dS_norm += vars.h(i,j) / vars.chi() * dS_U(i) * dS_U(j); }
    dS_norm = sqrt(dS_norm);
    FOR(i) { dS_U(i) /= dS_norm; }

    Tensor::Rank1 dS_L{};
    FOR(i)
    {
        // dS_L[i] = dS_U[i];
        dS_L(i)=0;
        FOR(j) { dS_L(i) += vars.h(i,j) / vars.chi() * dS_U(j); }
    }

    amrex::Real Madm{};
    FOR(i, j, k)
    {
        Madm += dS_L(i) / (16. * M_PI * m_G_Newton) /
        (vars.chi() * sqrt(vars.chi())) * h_UU(i,k) *
        (vars.h(j,j) * d1_chi(k) - vars.chi() * d1_h(j,j,k));
        FOR(l)
        {
            Madm += dS_L(i) / (16. * M_PI * m_G_Newton) /
                (vars.chi() * sqrt(vars.chi())) * h_UU(j,l) * h_UU(i,k) *
                (vars.chi() * d1_h(k,l,j) - vars.h(k,l) * d1_chi(j));
        }

    
    }
    // assign values of ADM Mass in output box
    const auto adm_cell_data = adm_state.cellData(ix, iy, iz);
    adm_cell_data[m_c_Madm] = Madm;
}


void ADMQuantities::set_up(int a_state_index)
{
    int num_ghosts  = 2; // no advection terms so only need 2 ghost cells

    const auto &comp_names = var_names;
    auto &derive_lst       = amrex::AmrLevel::get_derive_lst();
    const auto &desc_lst   = amrex::AmrLevel::get_desc_lst();

    derive_lst.add(
        name, amrex::IndexType::TheCellType(),
        static_cast<int>(comp_names.size()), comp_names, compute_mf,
        [=](const amrex::Box &box) { return amrex::grow(box, num_ghosts); },
        &amrex::cell_quartic_interp);

    // Get all the vars to allow us to use the CCZ4Vars class
    // Technically not all needed but probably doesn't hurt performance
    derive_lst.addComponent(name, desc_lst, a_state_index, 0, NUM_CCZ4_VARS);
}

void ADMQuantities::compute_mf(amrex::MultiFab &out_mf, int dcomp, int ncomp,
                             const amrex::MultiFab &src_mf,
                             const amrex::Geometry &geomdata,
                             amrex::Real /*time*/, const int * /*bcrec*/,
                             int /*level*/)
{
    const auto &out_arrays = out_mf.arrays();
    const auto &src_arrays = src_mf.const_arrays();
    int imadm               = dcomp;
    // Interval imom          = Interval(dcomp + 1, dcomp + AMREX_SPACEDIM);
    ADMQuantities admquantities(geomdata.CellSize(0),imadm);

    amrex::ParallelFor(
        out_mf, out_mf.nGrowVect(),
        [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz) noexcept
        { admquantities(ix, iy, iz, out_arrays[box_no], src_arrays[box_no]); });
}

#endif /* ADMQUANTITIES_IMPL_HPP_ */
