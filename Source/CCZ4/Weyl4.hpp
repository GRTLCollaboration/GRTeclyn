/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef WEYL4_HPP_
#define WEYL4_HPP_

#include "CCZ4Geometry.hpp"
#include "CCZ4RHS.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRParmParse.hpp"
#include "StateVariables.hpp" //This files needs c_NUM - total number of components
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include <array>

// AMReX Includes
#include <AMReX_AmrLevel.H>

//! Struct for the E and B fields
struct EBFields_t
{
    Tensor<2, amrex::Real> E; //!< Electric component of Weyltensor
    Tensor<2, amrex::Real> B; //!< Magnetic component of Weyltensor
};

//! Struct for the null tetrad
struct Tetrad_t
{
    Tensor<1, amrex::Real> u; //!< the vector u^i
    Tensor<1, amrex::Real> v; //!< the vector v^i
    Tensor<1, amrex::Real> w; //!< the vector w^i
};

//! Struct for the Newman Penrose scalar
struct weyl_scalar_t
{
    amrex::Real Real; // Real component
    amrex::Real Im;   // Imaginary component
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

    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    //! Constructor of class Weyl4
    /*!
        Takes in the centre for the calculation of the tetrads, grid spacing and
        the formulation.
    */
    // TODO: Remove dependence on formulation?
    Weyl4(const std::array<double, AMREX_SPACEDIM> &a_center, double a_dx,
          int a_out_comp, int a_formulation = CCZ4RHS<>::USE_CCZ4)
        : m_center(a_center), m_dx(a_dx), m_deriv(a_dx), m_out_comp(a_out_comp),
          m_formulation(a_formulation)
    {
    }
    // NOLINTEND(bugprone-easily-swappable-parameters)

    //! Computes Weyl4 in a cell
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    operator()(int ix, int iy, int iz,
               const amrex::Array4<amrex::Real> &weyl_scalars,
               const amrex::Array4<amrex::Real const> &state) const;

    AMREX_FORCE_INLINE static void set_up(int a_state_index);

    // Has signature of DeriveFuncMF so that it can be stored in the derive_lst
    AMREX_FORCE_INLINE static void
    compute_mf(amrex::MultiFab &out_mf, int out_comp, int ncomp,
               const amrex::MultiFab &src_mf, const amrex::Geometry &geomdata,
               amrex::Real /*time*/, const int * /*bcrec*/, int /*level*/);

  protected:
    std::array<double, AMREX_SPACEDIM> m_center; //!< The grid center
    double m_dx;                                 //!< the grid spacing
    FourthOrderDerivatives m_deriv; //!< for calculating derivs of vars
    int m_out_comp;    //!< Which commponent to store Weyl4_Re (Weyl4_Im will be
                       //!< m_out_comp+1)
    int m_formulation; //!< CCZ4 or BSSN?

    //! Compute spatial volume element
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor<3, amrex::Real>
    compute_epsilon3_LUU(const CCZ4Vars &vars,
                         const Tensor<2, amrex::Real> &h_UU) const;

    //! Calculation of Weyl_4 scalar
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE weyl_scalar_t
    compute_Weyl4(const EBFields_t &ebfields, const CCZ4Vars &vars,
                  const Tensor<2, amrex::Real> &h_UU,
                  const Coordinates &coords) const;

    //! Calculation of the tetrads
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tetrad_t
    compute_null_tetrad(const CCZ4Vars &vars,
                        const Tensor<2, amrex::Real> &h_UU,
                        const Coordinates &coords) const;

    //! Calulation of the decomposition of the Weyl tensor in Electric and
    //! Magnetic fields
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE EBFields_t
    compute_EB_fields(const CCZ4Vars &vars, const CCZ4D1Vars &d1,
                      const Tensor<2, amrex::Real> &d2_chi,
                      const Tensor<4, amrex::Real> &d2_h,
                      const Tensor<3, amrex::Real> &epsilon3_LUU,
                      const Tensor<2, amrex::Real> &h_UU,
                      const chris_t &chris) const;
};

#include "Weyl4.impl.hpp"

#endif /* WEYL4_HPP_ */
