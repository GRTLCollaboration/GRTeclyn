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

    //    AMREX_ENUM_IN_CLASS(Options, Mean, Variance);

    //! Constructor
    SpectralWishes(int state_index, const amrex::Real time, int lev, int ngrow)
        : m_time(time), m_lev(lev), m_state_index(state_index), m_ngrow(ngrow)
    {
        fill_params();
    };

    void fill_params(void)
    {
        GRParmParse spectral_wishes_pp("spectral_wishes");

        spectral_wishes_pp.query("operation", m_operation);
        spectral_wishes_pp.query("variable_name", m_var_name);
    };
    //! Main compute function wrapper (that calls specific compute function)
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real compute(GRAmr *gramr_ptr);

    //! Compute members

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    compute_mean(GRAmr *gramr_ptr);

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    compute_variance(GRAmr *gramr_ptr);

    //! Helper function for derived quantities
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real get_ncomp();

  protected:
    matter_t m_matter;        //!< The matter object, e.g. a scalar field
    const amrex::Real m_time; // Time as measured by AmrLevel
    const int m_lev;          // AMR level
    int m_state_index;
    const int
        m_ngrow; // Number of ghost cells to use in Spectral Wish calculation,
                 // should match what is used in the state MF in most cases.
    std::string m_var_name;  // The variable you wish to operate on
    std::string m_operation; // The operation e.g. mean, variance.
};

#include "SpectralWishes.impl.hpp"

#endif /* SPECTRALWISHES_HPP_ */
