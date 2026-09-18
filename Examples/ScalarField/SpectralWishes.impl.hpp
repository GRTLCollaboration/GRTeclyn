/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// Calculate the desired quantity

#ifndef SPECTRALWISHES_IMPL_HPP_
#define SPECTRALWISHES_IMPL_HPP_

#include "SpectralWishes.hpp"

template <class matter_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
SpectralWishes<matter_t>::compute_mean(const amrex::Geometry &geom,
                                       const amrex::MultiFab &src_mf, int ncomp)
{
    amrex::Real mean{0.};

    // Sum all the values in ncomp
    AMREX_ASSERT(ncomp < NUM_VARS);

    mean = src_mf.sum(ncomp);

    // Calculate total number of cells on this level
    // I think this is the same for all ranks
    const auto problo = geom.ProbLo();
    const auto probhi = geom.ProbHi();
    const auto dx     = geom.CellSizeArray();
    int n_cells       = AMREX_D_TERM((probhi[0] - problo[0]) / dx[0],
                                     +((probhi[1] - problo[1]) / dx[1]),
                                     +((probhi[2] - problo[2]) / dx[2]));

    mean /= static_cast<amrex::Real>(n_cells);

    amrex::AllPrint() << "Mean on Rank " << amrex::ParallelDescriptor::MyProc()
                      << ": " << mean << "\n";
    amrex::AllPrint() << "Number of cells on Rank "
                      << amrex::ParallelDescriptor::MyProc() << ": " << n_cells
                      << "\n";

    return mean;
}

// template <class matter_t>
// void SpectralWishes<matter_t>::set_up(int a_state_index, bool
// a_calc_mom_norm)
// {

//     int num_ghosts = 2;

//     int ncomp = 1;

//     auto &derive_lst     = amrex::AmrLevel::get_derive_lst();
//     const auto &desc_lst = amrex::AmrLevel::get_desc_lst();

//     // Add Constraints to the derive list
//     derive_lst.add(
//         SpectralWishes::class_name, amrex::IndexType::TheCellType(), ncomp,
//         SpectralWishes::class_name, SpectralWishes::compute_mf,
//         [=](const amrex::Box &box) { return amrex::grow(box, num_ghosts); },
//         &amrex::cell_quartic_interp);

//     derive_lst.addComponent(SpectralWishes::class_name, desc_lst,
//     a_state_index,
//                             0, NUM_VARS);
// }
// template <class matter_t>
// void SpectralWishes<matter_t>::compute_mf(amrex::MultiFab &out_mf, int dcomp,
//                                           int ncomp,
//                                           const amrex::MultiFab &src_mf,
//                                           const amrex::Geometry &geomdata,
//                                           amrex::Real /*time*/,
//                                           const int * /*bcrec*/, int
//                                           /*level*/)
// {
//     const auto &out_arrays = out_mf.arrays();
//     const auto &src_arrays = src_mf.const_arrays();

//     amrex::Real dx = geomdata.CellSize(0);

//     SpectralWishes<matter_t> mySpectralWishes;

//     amrex::ParallelFor(
//         out_mf, out_mf.nGrowVect(),
//         [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz) noexcept
//         { constraints(ix, iy, iz, out_arrays[box_no], src_arrays[box_no]);
//         });
// }

#endif /* SPECTRALWISHES_IMPL_HPP_ */
