/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(EMTENSOR_HPP)
#error "This file should only be included through EMTensor.hpp"
#endif

#ifndef EMTENSOR_IMPL_HPP
#define EMTENSOR_IMPL_HPP

#include <AMReX_AmrLevel.H>

#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Interval.hpp"
#include "simd.hpp"

template <class matter_t, enum EMTensorOptions em_tensor_options>
amrex::Vector<std::string> EMTensor<matter_t, em_tensor_options>::var_names()
{
    if constexpr (em_tensor_options == EMTensorOptions::justEnergyDensity)
    {
        return {"rho"};
    }
    else if (em_tensor_options == EMTensorOptions::energyAndMomentumDensities)
    {
        return {"rho", "j_1", "j_2", "j_3"};
    }
    else
    {
        return {"rho",  "j_1",  "j_2",  "j_3",  "S_11",
                "S_12", "S_13", "S_22", "S_23", "S_33"};
    }
}

template <class matter_t, enum EMTensorOptions em_tensor_options>
EMTensor<matter_t, em_tensor_options>::EMTensor(double a_dx, int a_dcomp)
    : m_deriv(a_dx), m_dcomp(a_dcomp)
{
}

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
template <class matter_t, enum EMTensorOptions em_tensor_options>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
EMTensor<matter_t, em_tensor_options>::compute(
    int i, int j, int k, const amrex::Array4<amrex::Real> &out_arrays,
    const amrex::Array4<const amrex::Real> &in_arrays) const
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    const auto vars = load_vars<Vars>(in_arrays.cellData(i, j, k));
    const auto d1   = m_deriv.template diff1<Vars>(i, j, k, in_arrays);

    using namespace TensorAlgebra;

    const auto h_UU  = compute_inverse_sym(vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);

    const auto emtensor = m_matter.compute_emtensor(vars, d1, h_UU, chris.ULL);

    out_arrays(i, j, k, m_dcomp) = emtensor.rho;

    if constexpr (em_tensor_options ==
                      EMTensorOptions::energyAndMomentumDensities ||
                  em_tensor_options == EMTensorOptions::allDensities)
    {
#if DEFAULT_TENSOR_DIM == 3
        FOR (i)
        {
            out_arrays(i, j, k, m_dcomp + 1 + i) = emtensor.j(i);
        }
#endif
    }

    if constexpr (em_tensor_options == EMTensorOptions::allDensities)
    {
#if DEFAULT_TENSOR_DIM == 3
        out_arrays(i, j, k, m_dcomp + 4) = emtensor.S(0, 0);
        out_arrays(i, j, k, m_dcomp + 5) = emtensor.S(0, 1);
        out_arrays(i, j, k, m_dcomp + 6) = emtensor.S(0, 2);
        out_arrays(i, j, k, m_dcomp + 7) = emtensor.S(1, 1);
        out_arrays(i, j, k, m_dcomp + 8) = emtensor.S(1, 2);
        out_arrays(i, j, k, m_dcomp + 9) = emtensor.S(2, 2);
    }
}
#endif

template <class matter_t, enum EMTensorOptions em_tensor_options>
AMREX_FORCE_INLINE void
EMTensor<matter_t, em_tensor_options>::set_up(int a_state_index)
{
    int num_ghosts = 2;

    auto derive_var_names = var_names();

    auto &derive_lst     = amrex::AmrLevel::get_derive_lst();
    const auto &desc_lst = amrex::AmrLevel::get_desc_lst();

    // Add EMTensor to the derive list
    derive_lst.add(
        name, amrex::IndexType::TheCellType(), derive_var_names.size(),
        derive_var_names, EMTensor::compute_mf, [=](const amrex::Box &box)
        { return amrex::grow(box, num_ghosts); }, &amrex::cell_quartic_interp);

    derive_lst.addComponent(name, desc_lst, a_state_index, 0, NUM_VARS);
}

template <class matter_t, enum EMTensorOptions em_tensor_options>
AMREX_FORCE_INLINE void EMTensor<matter_t, em_tensor_options>::compute_mf(
    amrex::MultiFab &out_mf, int dcomp, int ncomp,
    const amrex::MultiFab &src_mf, const amrex::Geometry &geomdata,
    amrex::Real /*time*/, const int * /*bcrec*/, int /*level*/)
{
    const auto &out_arrays = out_mf.arrays();
    const auto &src_arrays = src_mf.const_arrays();

    EMTensor<matter_t, em_tensor_options> em_tensor(geomdata.CellSize(0),
                                                    dcomp);

    amrex::ParallelFor(
        out_mf,
        [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
        {
            em_tensor.compute(i, j, k, out_arrays[box_no], src_arrays[box_no]);
        });
}

#endif /* EMTENSOR_IMPL_HPP */
