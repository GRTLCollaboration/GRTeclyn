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

// struct SpectralOptions
// {

//     AMREX_ENUM_IN_CLASS(Op, Mean, Variance);

//     static std::string get_name(Op op)
//     {
//         return amrex::getEnumNameString(hero);
//     }

//     static Op get_op(std::string_view name) { return
//     amrex::getEnum<Op>(name); }

//     Op m_operation;
// };

template <class matter_t> class SpectralWishes
{
  public:

    //    AMREX_ENUM_IN_CLASS(Options, Mean, Variance);

    //! Constructor
    SpectralWishes(int state_index, const amrex::Real time, int lev, int ngrow)
    {
        m_time        = time;
        m_lev         = lev;
        m_state_index = state_index;
        m_ngrow       = ngrow;
        fill_params();
    };

    void fill_params(void)
    {
        GRParmParse scalar_field_pp("scalar_field");

        scalar_field_pp.query("wish_op", m_operation);
        scalar_field_pp.query("wish_var", m_var_name);
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
    matter_t m_matter;  //!< The matter object, e.g. a scalar field
    amrex::Real m_time; // Time as measured by AmrLevel
    int m_lev;          // AMR level
    int m_state_index;
    int m_ngrow; // Number of ghost cells to use in Spectral Wish calculation,
                 // should match what is used in the state MF in most cases.
    std::string m_var_name;  // The variable you wish to operate on
    std::string m_operation; // The operation e.g. mean, variance.
};

#include "SpectralWishes.impl.hpp"

#endif /* SPECTRALWISHES_HPP_ */
