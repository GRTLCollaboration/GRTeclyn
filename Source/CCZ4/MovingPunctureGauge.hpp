/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef MOVINGPUNCTUREGAUGE_HPP_
#define MOVINGPUNCTUREGAUGE_HPP_

#include "CCZ4AdvecVars.hpp"
#include "CCZ4D1Vars.hpp"
#include "CCZ4D2Vars.hpp"
#include "CCZ4Vars.hpp"
#include "DimensionDefinitions.hpp"
#include <AMReX_Array.H>
#include <AMReX_GpuQualifiers.H>
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
    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    rhs_gauge(const amrex::CellData<amrex::Real> &rhs_cell_data,
              const CCZ4Vars &vars, const amrex::Real &advec_lapse,
              const TensorArray::Rank1 &advec_shift,
              const TensorArray::Rank1 &advec_B,
              const TensorArray::Rank1 &advec_Gamma) const
    // NOLINTEND(bugprone-easily-swappable-parameters)
    {
        rhs_cell_data[c_lapse] = m_params.lapse_advec_coeff * advec_lapse -
                                 m_params.lapse_coeff *
                                     pow(vars.lapse(), m_params.lapse_power) *
                                     (vars.K() - 2.0 * vars.Theta());

        FOR (i)
        {
            rhs_cell_data[c_shift1 + i] =
                m_params.shift_advec_coeff * advec_shift(i) +
                m_params.shift_Gamma_coeff * vars.B(i);
            rhs_cell_data[c_B1 + i] =
                m_params.shift_advec_coeff * advec_B(i) -
                m_params.shift_advec_coeff * advec_Gamma(i) +
                rhs_cell_data[c_Gamma1 + i] - m_params.eta * vars.B(i);
        }
    }
};

#endif /* MOVINGPUNCTUREGAUGE_HPP_ */
