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

#include "GRParmParse.hpp"

template <class matter_t> class SpectralWishes
{
  public:

    //! Constructor
    SpectralWishes(int state_index, const amrex::Real time, int lev)
    {
        m_time        = time;
        m_lev         = lev;
        m_state_index = state_index;
        fill_params();
    };

    void fill_params(void)
    {
        GRParmParse scalar_field_pp("scalar_field");

        scalar_field_pp.query("wish", m_var_name);
    };
    //! Compute members
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    compute_mean(GRAmr *gramr_ptr, const amrex::MultiFab &src_mf);

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    compute_variance(GRAmr *gramr_ptr, const amrex::MultiFab &src_mf);

    //! Helper function for derived quantities
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE std::unique_ptr<amrex::MultiFab>
    get_derived_mf(GRAmr *gramr_ptr, int &ncomp, const int &ngrow);

    // static void set_up(int a_state_index);

    // // Has signature of DeriveFuncMF so that it can be stored in the
    // derive_lst static void compute_mf(amrex::MultiFab &out_mf, int dcomp,
    // int ncomp,
    //                        const amrex::MultiFab &src_mf,
    //                        const amrex::Geometry &geomdata,
    //                        amrex::Real /*time*/, const int * /*bcrec*/,
    //                        int /*level*/);

    static inline const std::string class_name = "SpectralWishes";

  protected:
    matter_t m_matter;  //!< The matter object, e.g. a scalar field
    amrex::Real m_time; // Time as measured by AmrLevel
    int m_lev;          // AMR level
    int m_state_index;
    std::string m_var_name;  // The variable you wish to operate on
    std::string m_operation; // The operation e.g. mean, variance.
};

#include "SpectralWishes.impl.hpp"

#endif /* SPECTRALWISHES_HPP_ */
