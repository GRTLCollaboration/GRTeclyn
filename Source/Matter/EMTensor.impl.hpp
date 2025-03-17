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

template <class matter_t>
EMTensor<matter_t>::EMTensor(const double dx, const int a_c_rho,
                             const Interval a_c_Si, const Interval a_c_Sij)
    : m_deriv(dx), m_c_rho(a_c_rho), m_c_Si(a_c_Si), m_c_Sij(a_c_Sij)
{
    if (m_c_Si.size() != 0)
    {
        // Si is a vector
        AMREX_ASSERT(m_c_Si.size() == DEFAULT_TENSOR_DIM);
    }

    if (m_c_Sij.size() != 0)
    {
        // Sij is a symmetric tensor
        AMREX_ASSERT(m_c_Sij.size() ==
                     DEFAULT_TENSOR_DIM * (DEFAULT_TENSOR_DIM + 1) / 2);
    }
}

template <class matter_t>
template <class data_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
EMTensor<matter_t>::compute(int i, int j, int k,
                            const amrex::Array4<data_t> &out_mf,
                            const amrex::Array4<const data_t> &in_mf) const
{
    const auto vars = load_vars<Vars>(in_mf.cellData(i, j, k));
    const auto d1   = m_deriv.template diff1<Vars>(i, j, k, in_mf);

    using namespace TensorAlgebra;

    const auto h_UU  = compute_inverse_sym(vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);

    const auto emtensor = m_matter.compute_emtensor(vars, d1, h_UU, chris.ULL);

    if (m_c_rho >= 0)
    {
        out_mf(i, j, k, m_c_rho) = emtensor.rho;
    }

    if (m_c_Si.size() > 0)
    {
#if DEFAULT_TENSOR_DIM == 3
        FOR (i)
        {
            out_mf(i, j, k, m_c_Si.begin() + i) = emtensor.Si[i];
        }
#endif
    }

    if (m_c_Sij.size() > 0)
    {
#if DEFAULT_TENSOR_DIM == 3
        out_mf(i, j, k, m_c_Sij.begin())     = emtensor.Sij[0][0];
        out_mf(i, j, k, m_c_Sij.begin() + 1) = emtensor.Sij[0][1];
        out_mf(i, j, k, m_c_Sij.begin() + 2) = emtensor.Sij[0][2];
        out_mf(i, j, k, m_c_Sij.begin() + 3) = emtensor.Sij[1][1];
        out_mf(i, j, k, m_c_Sij.begin() + 4) = emtensor.Sij[1][2];
        out_mf(i, j, k, m_c_Sij.begin() + 5) = emtensor.Sij[2][2];

#endif
    }
}

template <class matter_t>
AMREX_FORCE_INLINE void EMTensor<matter_t>::set_up(int a_state_index,
                                                   bool do_all_components)
{

    if (do_all_components)
    {
        var_names = ArrayTools::concatenate(var_names, extra_var_names);
    }

    int num_ghosts = 2;

    auto &derive_lst     = amrex::AmrLevel::get_derive_lst();
    const auto &desc_lst = amrex::AmrLevel::get_desc_lst();

    // Add EMTensor to the derive list
    derive_lst.add(
        name, amrex::IndexType::TheCellType(), var_names.size(), var_names,
        EMTensor::compute_mf, [=](const amrex::Box &box)
        { return amrex::grow(box, num_ghosts); }, &amrex::cell_quartic_interp);

    derive_lst.addComponent(name, desc_lst, a_state_index, 0, NUM_VARS);
}

template <class matter_t>
AMREX_FORCE_INLINE void EMTensor<matter_t>::compute_mf(
    amrex::MultiFab &out_mf, int dcomp, int ncomp,
    const amrex::MultiFab &src_mf, const amrex::Geometry &geomdata,
    amrex::Real /*time*/, const int * /*bcrec*/, int /*level*/)
{
    const auto &out_arrays = out_mf.arrays();
    const auto &src_arrays = src_mf.const_arrays();

    // a_c_rho is stored starting from dcomp
    // a_c_Si is stored starting from rho (dcomp+1)
    // a_c_Sij is stored starting from a_c_Si

    int c_Si_begin{0}, c_Si_end{-1};
    int c_Sij_begin{0}, c_Sij_end{-1};

    // Depending on how many components are required (ncomp), also:

    // Compute the momentum density
    if (ncomp > 1)
    {
        c_Si_begin = dcomp + 1;
        c_Si_end   = c_Si_begin + DEFAULT_TENSOR_DIM;
    }

    // Compute the spatial stress-energy density
    if (ncomp > 3)
    {
        c_Sij_begin = c_Si_end;
        c_Sij_end =
            c_Sij_begin + DEFAULT_TENSOR_DIM * (DEFAULT_TENSOR_DIM + 1) / 2;
    }

    Interval my_c_Si(c_Si_begin, c_Si_end);
    Interval my_c_Sij(c_Sij_begin, c_Sij_end);

    EMTensor<matter_t> em_tensor(geomdata.CellSize(0), dcomp, my_c_Si,
                                 my_c_Sij);

    amrex::ParallelFor(
        out_mf,
        [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
        {
            em_tensor.compute(i, j, k, out_arrays[box_no], src_arrays[box_no]);
        });
}

#endif /* EMTENSOR_IMPL_HPP */
