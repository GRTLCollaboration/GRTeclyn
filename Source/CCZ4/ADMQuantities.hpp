/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// This compute class calculates ADM mass and angular momentum 

#ifndef ADMQUANTITIES_HPP_
#define ADMQUANTITIES_HPP_

// GRTeclyn includes
#include "CCZ4Geometry.hpp"
#include "CCZ4Vars.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRInterval.hpp"
#include "Interval.hpp"
#include "Coordinates.hpp"

#include "Tensor.hpp"
#include "TensorAlgebra.hpp"

// AMReX includes
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

// System includes
#include <array>

class ADMQuantities
{
  public:
    static inline const std::string name = "admquantities";

    /// Variable names
    static inline const amrex::Vector<std::string> var_names = {"Madm"};
    // static inline const amrex::Vector<std::string> var_names = {"Madm", "Jadm_1","Jadm_2", "Jadm_3"};

    // static inline const amrex::Vector<std::string> var_names_norm = {"ADM_mass", "ADM_mom"};
                                                                     

    /// Struct for ADM quantities
    struct adm_t
    {
        amrex::Real ADM_mass{};
        // Tensor::Rank1 mom{};
        // Tensor::Rank1 mom_abs_terms{};
    };

    // Constructor which allows specifying Ham and Mom vars
    // if the interval of a_c_Moms is of size 1, then
    // sqrt(Mom1^2 + Mom2^2 + Mom3^2) is stored in that variable
    // ...abs_terms stores the absolute value of the individual terms in the
    // conformally decomposed expressions which can be used in to normalize
    // the constraint violations
    // Any zero-length Interval or negative var is not calculated
    // If a positive interval is passed for one of a_c_Moms or
    // a_c_moms_abs_terms then it must have length consistent with
    // s_calc_mom_norm
    AMREX_FORCE_INLINE
    ADMQuantities(amrex::Real dx, int a_c_Madm);
    // , int a_c_Ham, const Interval &a_c_Moms,
    //             int a_c_Ham_abs_terms              = -1,
    //             const Interval &a_c_Moms_abs_terms = Interval(),
    //             amrex::Real cosmological_constant  = 0.0)

    AMREX_FORCE_INLINE AMREX_GPU_DEVICE void
    operator()(int ix, int iy, int iz,
               const amrex::Array4<amrex::Real> &adm_state,
               const amrex::Array4<amrex::Real const> &state) const;

    /// Adds the constraints to the derive list
    /// Call in variableSetUp()
    AMREX_FORCE_INLINE static void set_up(int a_state_index);

    // Has signature of DeriveFuncMF so that it can be stored in the derive_lst
    AMREX_FORCE_INLINE static void
    compute_mf(amrex::MultiFab &out_mf, int dcomp, int ncomp,
               const amrex::MultiFab &src_mf, const amrex::Geometry &geomdata,
               amrex::Real /*time*/, const int * /*bcrec*/, int /*level*/);

  protected:
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
    // static inline bool s_calc_mom_norm =
    //     false; // set to true with set_up() to store just sqrt(Mom1^2 + Mom2^2 +
    //            // Mom3^2) instead of Mom1, Mom2, Mom3 separately
    FourthOrderDerivatives m_deriv;
    // int m_c_Ham;
    int m_c_Madm;
    // Interval m_c_Moms;
    // int m_c_Ham_abs_terms = -1;
    // Interval m_c_Moms_abs_terms;
    std::array<amrex::Real, AMREX_SPACEDIM> m_center;
    amrex::Real m_G_Newton;
    amrex::Real m_dx;
};

#include "ADMQuantities.impl.hpp"

#endif /*ADMQUANTITIES_HPP_*/
