#ifndef RADIALRECIPE_MATTER_CONSTRAINTS_HPP_
#define RADIALRECIPE_MATTER_CONSTRAINTS_HPP_

#include "ConstraintsWithMatter.hpp"
#include "Interval.hpp"

#include <AMReX_MultiFab.H>

#include <array>

template <class matter_t>
void fill_matter_constraints(amrex::MultiFab &cst,
                             const amrex::MultiFab &state_new,
                             const matter_t &matter, amrex::Real dx0,
                             const std::array<double, AMREX_SPACEDIM> &center,
                             amrex::Real time)
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

#endif /* RADIALRECIPE_MATTER_CONSTRAINTS_HPP_ */
