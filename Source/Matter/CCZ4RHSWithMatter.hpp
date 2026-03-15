/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef CCZ4RHSWITHMATTER_HPP_
#define CCZ4RHSWITHMATTER_HPP_

#include "CCZ4Geometry.hpp"
#include "CCZ4RHS.hpp"
#include "CCZ4Vars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "MovingPunctureGaugeWithMatter.hpp"
#include "StateVariables.hpp" //This files needs NUM_VARS - total number of components
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"

//!  Calculates RHS using CCZ4 including matter terms, and matter variable
//!  evolution
/*!
     The class calculates the RHS evolution for all the variables. It inherits
   from the CCZ4RHS class, which it uses to do the non matter evolution of
   variables. It then adds in the additional matter terms to the CCZ4 evolution
   (those including the stress energy tensor), and calculates the evolution of
   the matter variables. It does not assume a specific form of matter but is
   templated over a matter class matter_t. Please see the class ScalarField as
   an example of a matter_t. \sa CCZ4RHS(), ScalarField()
*/

template <class matter_t, class gauge_t = MovingPunctureGaugeWithMatter,
          class deriv_t = FourthOrderDerivatives>
class CCZ4RHSWithMatter : public CCZ4RHS<gauge_t, deriv_t>
{
  public:
    // Use this alias for the same template instantiation as this class
    using CCZ4 = CCZ4RHS<gauge_t, deriv_t>;

    using params_t = CCZ4_params_t<typename gauge_t::params_t>;

    //!  Constructor of class MatterCCZ4
    /*!
       Inputs are the grid spacing, plus the CCZ4 evolution parameters and a
       matter object. It also takes the dissipation parameter sigma, and allows
       the formulation to be toggled between CCZ4 and BSSN. The default is CCZ4.
       It allows the user to set the value of Newton's constant, which is set to
       one by default.
    */
    CCZ4RHSWithMatter(params_t a_params, double a_dx, double a_sigma,
                      int a_formulation = CCZ4RHS<>::USE_CCZ4,
                      double a_G_Newton = 1.0,
                      std::array<double, AMREX_SPACEDIM> a_center = {0.0, 0.0, 0.0},
                      amrex::Real a_time = 0.0);

    CCZ4RHSWithMatter(matter_t a_matter, params_t a_params, double a_dx, double a_sigma,
                      int a_formulation = CCZ4RHS<>::USE_CCZ4,
                      double a_G_Newton = 1.0,
                      std::array<double, AMREX_SPACEDIM> a_center = {0.0, 0.0, 0.0},
                      amrex::Real a_time = 0.0);

    //!  The compute member which calculates the RHS at each point in the box
    //!  \sa matter_rhs_equation()
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    operator()(int ix, int iy, int iz,
               const amrex::Array4<amrex::Real> &rhs_state,
               const amrex::Array4<amrex::Real const> &state) const;

  protected:
    //! The function which adds in the EM Tensor terms to the CCZ4 rhs \sa
    //! compute()
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void add_emtensor_rhs(
        const amrex::CellData<amrex::Real>
            &rhs, //!< the RHS data for each variable at that point.
        const typename matter_t::Vars
            &state, //!< the value of the variables at the point.
        const typename matter_t::D1Vars
            &d1, //!< the value of the first derivatives of the variables.
        const Coordinates &coords
    ) const;

    // Class members
    matter_t m_matter; //!< The matter object, e.g. a scalar field.
    double m_G_Newton; //!< Newton's constant, set to one by default.
    double m_dx;
    std::array<double, AMREX_SPACEDIM> m_center;
    amrex::Real m_time;
};

#include "CCZ4RHSWithMatter.impl.hpp"

#endif /* CCZ4RHSWITHMATTER_HPP_ */
