/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef MOVINGPUNCTUREGAUGE_HPP_
#define MOVINGPUNCTUREGAUGE_HPP_

#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"

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

    template <class data_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    rhs_gauge(const amrex::CellData<data_t> &rhs_cell_data,
	      const data_t lapse,
	      const data_t K,
	      const data_t Theta,
	      const Tensor<1, data_t> B,
        const amrex::GpuArray<Tensor<1, data_t>, NUM_CCZ4_VARS>  & /*d1*/,
              //const diff2_vars_t<Tensor<2, data_t>> & /*d2*/,
              const amrex::GpuArray<data_t, NUM_CCZ4_VARS> &advec) const
    {
        rhs_cell_data[c_lapse] = m_params.lapse_advec_coeff * advec[c_lapse] -
                    m_params.lapse_coeff *
                        pow(lapse, m_params.lapse_power) *
                        (K - 2 * Theta);
        FOR (i)
        {
          rhs_cell_data[c_shift1 + i] = m_params.shift_advec_coeff * advec[c_shift1 + i] +
                           m_params.shift_Gamma_coeff * B[i];
          rhs_cell_data[c_B1 + i] = m_params.shift_advec_coeff * advec[c_B1 + i] -
                       m_params.shift_advec_coeff * advec[c_Gamma1 + i] +
                       rhs_cell_data[c_Gamma1 + i] - m_params.eta * B[i];
        }
    }
};

#endif /* MOVINGPUNCTUREGAUGE_HPP_ */
