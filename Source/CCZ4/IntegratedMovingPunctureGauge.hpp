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
class IntegratedMovingPunctureGauge : public MovingPunctureGauge<deriv_t>
{
  public:
    using base_t   = MovingPunctureGauge<deriv_t>;
    using params_t = typename base_t::params_t;

    IntegratedMovingPunctureGauge(double a_dx) : base_t(a_dx) {}

    /// Store the initial integrated Gamma-driver RHS in B.
    /** This makes the non-advective part of the initial shift RHS vanish. The
     * B field is subsequently frozen by calculate_gauge_rhs(), preserving the
     * subtraction throughout the evolution.
     */
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    operator()(int ix, int iy, int iz,
               const amrex::Array4<amrex::Real> &state) const
    {
        const amrex::CellData<amrex::Real> &state_cell_data =
            state.cellData(ix, iy, iz);
        const amrex::CellData<const amrex::Real> &const_state_cell_data =
            state_cell_data;
        const CCZ4Vars vars(const_state_cell_data);

        amrex::Real eta_of_x;
        this->compute_eta(eta_of_x, ix, iy, iz);

        FOR (i)
        {
            state_cell_data[c_B1 + i] =
                this->m_params.shift_Gamma_coeff * vars.Gamma(i) -
                eta_of_x * vars.shift(i);
        }
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

        const amrex::Real advec_lapse = this->m_deriv.advec_scalar(
            ix, iy, iz, state, shift_vector, c_lapse);
        const Tensor::Rank1 advec_shift = this->m_deriv.advec_vector(
            ix, iy, iz, state, shift_vector, c_shift1);

        amrex::Real eta_of_x;
        this->compute_eta(eta_of_x, ix, iy, iz);

        rhs_cell_data[c_lapse] =
            this->m_params.lapse_advec_coeff * advec_lapse -
            this->m_params.lapse_coeff *
                pow(vars.lapse(), this->m_params.lapse_power) *
                (vars.K() - 2.0 * vars.Theta());

        FOR (i)
        {
            rhs_cell_data[c_shift1 + i] =
                this->m_params.shift_advec_coeff * advec_shift(i) +
                this->m_params.shift_Gamma_coeff * vars.Gamma(i) -
                eta_of_x * vars.shift(i) - vars.B(i);
            rhs_cell_data[c_B1 + i] = 0.0;
        }
    }
};

#endif /* INTEGRATEDMOVINGPUNCTUREGAUGE_HPP_ */
