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
#include <AMReX_Array.H>

// This class only works for 3+1D
static_assert(GR_SPACEDIM == 3, "GR_SPACEDIM must be 3");

//! Struct for the E and B fields
struct EBFields_t
{
    TensorArray::Rank2 E; //!< Electric component of Weyltensor
    TensorArray::Rank2 B; //!< Magnetic component of Weyltensor
};

//! Struct for the null tetrad
struct Tetrad_t
{
    TensorArray::Rank1 u; //!< the vector u^i
    TensorArray::Rank1 v; //!< the vector v^i
    TensorArray::Rank1 w; //!< the vector w^i
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
    compute_mf(amrex::MultiFab &out_mf, int dcomp, int ncomp,
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

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE TensorArray::Rank3
    compute_epsilon3_LUU(const CCZ4Vars &vars,
                         const TensorArray::Rank2 &h_UU) const;

    //! Calculation of Weyl_4 scalar
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE weyl_scalar_t
    compute_Weyl4(const EBFields_t &ebfields, const CCZ4Vars &vars,
                  const TensorArray::Rank2 &h_UU,
                  const Coordinates &coords) const;

    //! Calculation of the tetrads
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tetrad_t
    compute_null_tetrad(const CCZ4Vars &vars, const TensorArray::Rank2 &h_UU,
                        const Coordinates &coords) const;

    //! Calulation of the decomposition of the Weyl tensor in Electric and
    //! Magnetic fields
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE EBFields_t
    compute_EB_fields(const CCZ4Vars &vars, const CCZ4D1Vars &d1,
                      const TensorArray::Rank1Sym &d2_chi,
                      const TensorArray::Rank2Sym &d2_h,
                      const TensorArray::Rank3 &epsilon3_LUU,
                      const TensorArray::Rank2 &h_UU,
                      const chris_t &chris) const;
};

#include "Weyl4.impl.hpp"

#endif /* WEYL4_HPP_ */
