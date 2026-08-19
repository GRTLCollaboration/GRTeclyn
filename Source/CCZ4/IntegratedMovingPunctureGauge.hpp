/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef INTEGRATEDMOVINGPUNCTUREGAUGE_HPP_
#define INTEGRATEDMOVINGPUNCTUREGAUGE_HPP_

#include "CCZ4Vars.hpp"
#include "DimensionDefinitions.hpp"
#include "MovingPunctureGauge.hpp"

/// This is an example of a gauge class that can be used in the CCZ4RHS compute
/// class
/**
 * This class implements a slightly more generic version of the moving puncture
 * gauge. In particular it uses a Bona-Masso slicing condition of the form
 * f(lapse) = -c*lapse^(p-2)
 * and an Integrated version of the Gamma-driver shift condition
 * (see details in arXiv:gr-qc/0605030)
 **/
template <class deriv_t = FourthOrderDerivatives>
class IntegratedMovingPunctureGauge
{
  public:
    using params_t = typename MovingPunctureGauge<deriv_t>::params_t;

  protected:
    params_t m_params;
    deriv_t m_deriv;

  public:
    IntegratedMovingPunctureGauge(double a_dx) : m_deriv(a_dx)
    {
        m_params.fill_params();
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    calculate_gauge_rhs(int ix, int iy, int iz,
                        const amrex::Array4<amrex::Real> &rhs,
                        const amrex::Array4<const amrex::Real> &state) const
    {
        const amrex::CellData<amrex::Real> &rhs_cell_data =
            rhs.cellData(ix, iy, iz);
        const amrex::CellData<const amrex::Real> &state_cell_data =
            state.cellData(ix, iy, iz);
        const CCZ4Vars vars(state_cell_data);

        const Tensor::Rank1 shift_vector(
            {vars.shift(0), vars.shift(1), vars.shift(2)});

        const amrex::Real advec_lapse =
            m_deriv.advec_scalar(ix, iy, iz, state, shift_vector, c_lapse);
        const Tensor::Rank1 advec_shift =
            m_deriv.advec_vector(ix, iy, iz, state, shift_vector, c_shift1);

        rhs_cell_data[c_lapse] = m_params.lapse_advec_coeff * advec_lapse -
                                 m_params.lapse_coeff *
                                     pow(vars.lapse(), m_params.lapse_power) *
                                     (vars.K() - 2.0 * vars.Theta());

        FOR (i)
        {
            rhs_cell_data[c_shift1 + i] =
                m_params.shift_advec_coeff * advec_shift(i) +
                m_params.shift_Gamma_coeff * vars.Gamma(i) -
                m_params.eta * vars.shift(i) - vars.B(i);
            rhs_cell_data[c_B1 + i] = 0.0;
        }
    }
};

#endif /* INTEGRATEDMOVINGPUNCTUREGAUGE_HPP_ */
