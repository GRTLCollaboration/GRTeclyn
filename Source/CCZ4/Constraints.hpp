/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// This compute class calculates Hamiltonian and Momentum constraints

#ifndef CONSTRAINTS_HPP_
#define CONSTRAINTS_HPP_

// GRTeclyn includes
#include "CCZ4D1Vars.hpp"
#include "CCZ4Geometry2.hpp"
#include "Cell.hpp"
#include "ConstCCZ4Vars.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives2.hpp"
#include "GRInterval.hpp"
#include "Interval.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"

// AMReX includes
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

// System includes
#include <array>

class Constraints
{
  public:
    static inline const std::string name = "constraints";

    /// Variable names
    static inline const amrex::Vector<std::string> var_names = {"Ham", "Mom1",
                                                                "Mom2", "Mom3"};

    static inline const amrex::Vector<std::string> var_names_norm = {"Ham",
                                                                     "Mom"};

    /// Vars object for Constraints
    struct Vars
    {
        amrex::Real Ham{};
        amrex::Real Ham_abs_terms{};
        Tensor<1, amrex::Real> Mom;
        Tensor<1, amrex::Real> Mom_abs_terms;
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
    Constraints(double dx, int a_c_Ham, const Interval &a_c_Moms,
                int a_c_Ham_abs_terms              = -1,
                const Interval &a_c_Moms_abs_terms = Interval(),
                double cosmological_constant       = 0.0);

    AMREX_FORCE_INLINE AMREX_GPU_DEVICE void
    operator()(int ix, int iy, int iz,
               const amrex::Array4<amrex::Real> &constraints,
               const amrex::Array4<amrex::Real const> &state) const;

    /// Adds the constraints to the derive list
    /// Call in variableSetUp()
    AMREX_FORCE_INLINE static void set_up(int a_state_index,
                                          bool a_calc_mom_norm = false);

    // Has signature of DeriveFuncMF so that it can be stored in the derive_lst
    AMREX_FORCE_INLINE static void
    compute_mf(amrex::MultiFab &out_mf, int dcomp, int ncomp,
               const amrex::MultiFab &src_mf, const amrex::Geometry &geomdata,
               amrex::Real /*time*/, const int * /*bcrec*/, int /*level*/);

  protected:
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
    static inline bool s_calc_mom_norm =
        false; // set to true with set_up() to store just sqrt(Mom1^2 + Mom2^2 +
               // Mom3^2) instead of Mom1, Mom2, Mom3 separately
    FourthOrderDerivatives2 m_deriv;
    int m_c_Ham;
    Interval m_c_Moms;
    int m_c_Ham_abs_terms = -1;
    Interval m_c_Moms_abs_terms;
    double m_cosmological_constant;

    AMREX_GPU_DEVICE Vars constraint_equations(
        const ConstCCZ4Vars &vars, const CCZ4D1Vars &d1,
        const Tensor<2, amrex::Real> &d2_chi,
        const Tensor<4, amrex::Real> &d2_h, const Tensor<2, amrex::Real> &h_UU,
        const chris_t &chris) const;

    AMREX_FORCE_INLINE AMREX_GPU_DEVICE void
    store_vars(const Vars &out,
               const amrex::CellData<amrex::Real> &current_cell) const;
};

#include "Constraints.impl.hpp"

#endif /* CONSTRAINTS_HPP_ */
