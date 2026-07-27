#ifndef RADIALRECIPE_MATTER_CONSTRAINTS_HPP_
#define RADIALRECIPE_MATTER_CONSTRAINTS_HPP_

#include "ConstraintsWithMatter.hpp"
#include "Interval.hpp"

#include <AMReX_MultiFab.H>

#include <array>

//! Fill Hamiltonian / momentum constraints into comps 0-3.
//! The Constraints class forbids requesting Ham and Ham_abs in one pass, so
//! absolute-term fills go through ``fill_matter_constraint_abs_terms``.
template <class matter_t>
void fill_matter_constraints(amrex::MultiFab &cst,
                             const amrex::MultiFab &state_new,
                             const matter_t &matter, amrex::Real dx0,
                             const std::array<double, AMREX_SPACEDIM> &center,
                             amrex::Real time, bool /*with_abs_terms*/ = false)
{
    ConstraintsWithMatter<matter_t> my_constraints(
        matter, dx0, 1.0, 0, Interval(1, 3), center, time);
    for (amrex::MFIter mfi(cst, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.validbox();
        const auto arr       = cst.array(mfi);
        const auto src_arr   = state_new.const_array(mfi);
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int ix, int iy, int iz) noexcept
            { my_constraints(ix, iy, iz, arr, src_arr); });
    }
}

//! Second pass: write Ham_abs / Mom_abs into comps 4 / 5-7 of ``cst``.
//! ``cst`` must have at least 8 components. Constraints ctor requires Ham xor
//! Ham_abs, so this pass uses c_Ham=-1.
template <class matter_t>
void fill_matter_constraint_abs_terms(
    amrex::MultiFab &cst, const amrex::MultiFab &state_new,
    const matter_t &matter, amrex::Real dx0,
    const std::array<double, AMREX_SPACEDIM> &center, amrex::Real time)
{
    ConstraintsWithMatter<matter_t> my_constraints(
        matter, dx0, 1.0, /*c_Ham=*/-1, Interval(), center, time,
        /*c_Ham_abs=*/4, Interval(5, 7));
    for (amrex::MFIter mfi(cst, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.validbox();
        const auto arr       = cst.array(mfi);
        const auto src_arr   = state_new.const_array(mfi);
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int ix, int iy, int iz) noexcept
            { my_constraints(ix, iy, iz, arr, src_arr); });
    }
}

#endif /* RADIALRECIPE_MATTER_CONSTRAINTS_HPP_ */
