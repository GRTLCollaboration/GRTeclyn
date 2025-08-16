/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef MOVINGPUNCTUREGAUGE_HPP_
#define MOVINGPUNCTUREGAUGE_HPP_

#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include <AMReX_REAL.H>

/// This is an example of a gauge class that can be used in the CCZ4RHS compute
/// class
/**
 * This class implements a slightly more generic version of the moving puncture
 * gauge. In particular it uses a Bona-Masso slicing condition of the form
 * f(lapse) = -c*lapse^(p-2)
 * and a Gamma-driver shift condition
 **/
class MovingPunctureGauge
{
  public:
    struct params_t
    {
        // lapse params:
        double lapse_advec_coeff = 0.; //!< Switches advection terms in
                                       //! the lapse condition on/off
        double lapse_power = 1.; //!< The power p in \f$\partial_t \alpha = - c
                                 //!\alpha^p(K-2\Theta)\f$
        double lapse_coeff = 2.; //!< The coefficient c in \f$\partial_t \alpha
                                 //!= -c \alpha^p(K-2\Theta)\f$
        // shift params:
        double shift_Gamma_coeff = 0.75; //!< Gives the F in \f$\partial_t
                                         //!  \beta^i =  F B^i\f$
        double shift_advec_coeff = 0.;   //!< Switches advection terms in the
                                         //! shift condition on/off
        double eta = 1.; //!< The eta in \f$\partial_t B^i = \partial_t \tilde
                         //!\Gamma - \eta B^i\f$
    };

  protected:
    params_t m_params;

  public:
    MovingPunctureGauge(const params_t &a_params) : m_params(a_params) {}
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    rhs_gauge(CCZ4Vars &rhs, const ConstCCZ4Vars &vars, const CCZ4D1Vars &d1,
              const CCZ4D2Vars &d2, const CCZ4AdvecVars &advec) const
    {
        amrex::Real rhs_lapse;
        rhs_lapse = m_params.lapse_advec_coeff * advec.lapse -
                    m_params.lapse_coeff *
                        pow(vars.lapse(), m_params.lapse_power) *
                        (vars.K() - 2.0 * vars.Theta());
        rhs.store_lapse(rhs_lapse);

        Tensor<1, amrex::Real> rhs_shift;
        Tensor<1, amrex::Real> rhs_B;
        FOR (i)
        {
            rhs_shift[i] = m_params.shift_advec_coeff * advec.shift[i] +
                           m_params.shift_Gamma_coeff * vars.B(i);
            rhs_B[i] = m_params.shift_advec_coeff * advec.B[i] -
                       m_params.shift_advec_coeff * advec.Gamma[i] +
                       rhs.Gamma(i) - m_params.eta * vars.B(i);
        }
        rhs.store_shift(rhs_shift);
        rhs.store_B(rhs_B);
    }
};

#endif /* MOVINGPUNCTUREGAUGE_HPP_ */
