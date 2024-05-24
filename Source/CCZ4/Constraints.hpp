/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// This compute class calculates Hamiltonian and Momentum constraints

#ifndef CONSTRAINTS_HPP_
#define CONSTRAINTS_HPP_

#include "BSSNVars.hpp"
#include "Cell.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Interval.hpp"
#include "Tensor.hpp"
#include "simd.hpp"

#include "CCZ4Geometry.hpp"

#include <array>

class Constraints
{
  public:
    /// Variable names
    static inline const amrex::Vector<std::string> var_names = {"Ham", "Mom1",
                                                                "Mom2", "Mom3"};

    static inline const amrex::Vector<std::string> var_names_norm = {"Ham",
                                                                     "Mom"};

    /// CCZ4 variables
    template <class data_t> using MetricVars = BSSNVars::VarsNoGauge<data_t>;

    /// CCZ4 variables
    template <class data_t>
    using Diff2Vars = BSSNVars::Diff2VarsNoGauge<data_t>;

    /// Vars object for Constraints
    template <class data_t> struct Vars
    {
        data_t Ham{};
        data_t Ham_abs_terms{};
        Tensor<1, data_t> Mom;
        Tensor<1, data_t> Mom_abs_terms;
    };

    // Constructor which allows specifying Ham and Mom vars
    // if the interval of a_c_Moms is of size 1, then
    // sqrt(Mom1^2 + Mom2^2 + Mom3^2) is stored in that variable
    // ...abs_terms stores the absolute value of the individual terms in the
    // conformally decomposed expressions which can be used in to normalize
    // the constraint violations
    // Any zero-length Interval or negative var is not calculated
    Constraints(double dx, int a_c_Ham, const Interval &a_c_Moms,
                int a_c_Ham_abs_terms              = -1,
                const Interval &a_c_Moms_abs_terms = Interval(),
                double cosmological_constant       = 0.0);

    template <class data_t>
    AMREX_GPU_DEVICE void
    compute(int i, int j, int k, const amrex::Array4<data_t> &cst,
            const amrex::Array4<data_t const> &state) const;

  protected:
    FourthOrderDerivatives m_deriv;
    int m_c_Ham;
    Interval m_c_Moms;
    int m_c_Ham_abs_terms = -1;
    Interval m_c_Moms_abs_terms;
    double m_cosmological_constant;

    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t>
    AMREX_GPU_DEVICE Vars<data_t> constraint_equations(
        const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
        const diff2_vars_t<Tensor<2, data_t>> &d2,
        const Tensor<2, data_t> &h_UU, const chris_t<data_t> &chris) const;

    template <class data_t>
    AMREX_GPU_DEVICE void
    store_vars(const Vars<data_t> &out,
               const amrex::CellData<data_t> &current_cell) const;
};

#include "Constraints.impl.hpp"

#endif /* CONSTRAINTS_HPP_ */
