/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SPECTRALWISHES_HPP_
#define SPECTRALWISHES_HPP_

#include <AMReX_MultiFab.H>
#include <AMReX_ParallelContext.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>

template <class matter_t> class SpectralWishes
{
  public:

    //! Default Constructor for now
    SpectralWishes() = default;

    //! The compute member which calculates the constraints at each point in the
    //! box
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    compute_mean(const amrex::Geometry &geom, const amrex::MultiFab &src_mf,
                 int ncomp);

    // static void set_up(int a_state_index);

    // // Has signature of DeriveFuncMF so that it can be stored in the
    // derive_lst static void compute_mf(amrex::MultiFab &out_mf, int dcomp, int
    // ncomp,
    //                        const amrex::MultiFab &src_mf,
    //                        const amrex::Geometry &geomdata,
    //                        amrex::Real /*time*/, const int * /*bcrec*/,
    //                        int /*level*/);

    static inline const std::string class_name = "SpectralWishes";

  protected:
    matter_t m_matter; //!< The matter object, e.g. a scalar field
};

#include "SpectralWishes.impl.hpp"

#endif /* SPECTRALWISHES_HPP_ */
