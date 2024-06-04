/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef WEYL4_HPP_
#define WEYL4_HPP_

#include "CCZ4RHS.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "StateVariables.hpp" //This files needs c_NUM - total number of components
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "simd.hpp"
#include <array>

//! Struct for the E and B fields
template <class data_t> struct EBFields_t
{
    Tensor<2, data_t> E; //!< Electric component of Weyltensor
    Tensor<2, data_t> B; //!< Magnetic component of Weyltensor
};

//! Struct for the null tetrad
template <class data_t> struct Tetrad_t
{
    Tensor<1, data_t> u; //!< the vector u^i
    Tensor<1, data_t> v; //!< the vector v^i
    Tensor<1, data_t> w; //!< the vector w^i
};

//! Struct for the Newman Penrose scalar
template <class data_t> struct NPScalar_t
{
    data_t Real; // Real component
    data_t Im;   // Imaginary component
};

//!  Calculates the Weyl4 scalar for spacetimes without matter content
/*!
   This class calculates the Weyl4 scalar real and im parts using definitions
   from Alcubierres book "Introduction to 3+1 Numerical Relativity". We use a
   decomposition of the Weyl tensor in electric and magnetic parts, which then
   is used together with the tetrads defined in "gr-qc/0104063" to calculate the
   Weyl4 scalar.
*/
class Weyl4
{
  public:
    /// derive record name
    static inline const std::string name = "Weyl4";

    /// Variable names
    static inline const amrex::Vector<std::string> var_names = {"Weyl4_Re",
                                                                "Weyl4_Im"};

    // Use the variable definitions containing the needed quantities
    template <class data_t> using Vars = CCZ4Vars::VarsWithGauge<data_t>;
    template <class data_t>
    using Diff2Vars = ADMConformalVars::Diff2VarsNoGauge<data_t>;

    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    //! Constructor of class Weyl4
    /*!
        Takes in the centre for the calculation of the tetrads, grid spacing and
        the formulation.
    */
    Weyl4(const std::array<double, AMREX_SPACEDIM> &a_center, double a_dx,
          int a_dcomp, int a_formulation = CCZ4RHS<>::USE_CCZ4)
        : m_center(a_center), m_dx(a_dx), m_deriv(a_dx), m_dcomp(a_dcomp),
          m_formulation(a_formulation)
    {
    }
    // NOLINTEND(bugprone-easily-swappable-parameters)

    //! Computes Weyl4 in a cell
    template <class data_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    compute(int i, int j, int k, const amrex::Array4<data_t> &a_derive_array,
            const amrex::Array4<data_t const> &a_state_array) const;

    static void set_up(int a_state_index);

  protected:
    std::array<double, AMREX_SPACEDIM> m_center; //!< The grid center
    double m_dx;                                 //!< the grid spacing
    FourthOrderDerivatives m_deriv; //!< for calculating derivs of vars
    int m_dcomp;       //!< Which commponent to store Weyl4_Re (Weyl4_Im will be
                       //!< m_dcomp+1)
    int m_formulation; //!< CCZ4 or BSSN?

    //! Compute spatial volume element
    template <class data_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor<3, data_t>
    compute_epsilon3_LUU(const Vars<data_t> &vars,
                         const Tensor<2, data_t> &h_UU) const;

    //! Calculation of Weyl_4 scalar
    template <class data_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE NPScalar_t<data_t>
    compute_Weyl4(const EBFields_t<data_t> &ebfields, const Vars<data_t> &vars,
                  const Vars<Tensor<1, data_t>> &d1,
                  const Diff2Vars<Tensor<2, data_t>> &d2,
                  const Tensor<2, data_t> &h_UU,
                  const Coordinates<data_t> &coords) const;

    //! Calculation of the tetrads
    template <class data_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tetrad_t<data_t>
    compute_null_tetrad(const Vars<data_t> &vars, const Tensor<2, data_t> &h_UU,
                        const Coordinates<data_t> &coords) const;

    //! Calulation of the decomposition of the Weyl tensor in Electric and
    //! Magnetic fields
    template <class data_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE EBFields_t<data_t> compute_EB_fields(
        const Vars<data_t> &vars, const Vars<Tensor<1, data_t>> &d1,
        const Diff2Vars<Tensor<2, data_t>> &d2,
        const Tensor<3, data_t> &epsilon3_LUU, const Tensor<2, data_t> &h_UU,
        const chris_t<data_t> &chris) const;
};

#include "Weyl4.impl.hpp"

#endif /* WEYL4_HPP_ */
