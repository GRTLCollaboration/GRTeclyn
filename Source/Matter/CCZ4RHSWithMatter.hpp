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
#include "FourthOrderDerivatives.hpp"
#include "MovingPunctureGaugeWithMatter.hpp"
#include "StateVariables.hpp" //This files needs NUM_VARS - total number of components
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
                      double a_G_Newton = 1.0);

    CCZ4RHSWithMatter(matter_t a_matter, params_t a_params, double a_dx,
                      double a_sigma, int a_formulation = CCZ4RHS<>::USE_CCZ4,
                      double a_G_Newton = 1.0);

    //! Add the stress-energy tensor terms to the CCZ4 and gauge RHS.
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void add_emtensor_rhs(
        const int ix, const int iy, const int iz,
        const amrex::Array4<amrex::Real>
            &rhs_state, //!< the RHS data for each variable at that point.
        const amrex::Array4<const amrex::Real>
            &state) //!< the current value of the variables at the point.
        const;

    //! Add the evolution equations for the matter variables.
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    add_matter_rhs(int ix, int iy, int iz,
                   const amrex::Array4<amrex::Real> &rhs_state,
                   const amrex::Array4<const amrex::Real> &state) const;

    //! Add dissipation to the CCZ4 and matter variables.
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    apply_dissipation(int ix, int iy, int iz,
                      const amrex::Array4<amrex::Real> &rhs_state,
                      const amrex::Array4<const amrex::Real> &state) const;

  protected:
    // Class members
    matter_t m_matter; //!< The matter object, e.g. a scalar field.
    double m_G_Newton; //!< Newton's constant, set to one by default.
};

#include "CCZ4RHSWithMatter.impl.hpp"

#endif /* CCZ4RHSWITHMATTER_HPP_ */
