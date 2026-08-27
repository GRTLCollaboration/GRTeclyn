/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef MOVINGPUNCTUREGAUGE_HPP_
#define MOVINGPUNCTUREGAUGE_HPP_

#include "CCZ4Vars.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRParmParse.hpp"
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
template <class deriv_t = FourthOrderDerivatives> class MovingPunctureGauge
{
  public:
    struct params_t
    {
        // lapse params:
        amrex::Real lapse_advec_coeff; //!< Switches advection terms in
                                       //! the lapse condition on/off
        amrex::Real lapse_power; //!< The power p in \f$\partial_t \alpha = - c
                                 //!\alpha^p(K-2\Theta)\f$
        amrex::Real lapse_coeff; //!< The coefficient c in \f$\partial_t \alpha
                                 //!= -c \alpha^p(K-2\Theta)\f$
        // shift params:
        amrex::Real shift_Gamma_coeff; //!< Gives the F in \f$\partial_t
                                       //!  \beta^i =  F B^i\f$
        amrex::Real shift_advec_coeff; //!< Switches advection terms in the
                                       //! shift condition on/off
        amrex::Real eta; //!< The eta in \f$\partial_t B^i = \partial_t \tilde
                         //!\Gamma - \eta B^i\f$

        static void check_params()
        {
            GRParmParse gauge_pp("gauge");
            // Lapse evolution
            amrex::Real lapse_advec_coeff = 1.;
            gauge_pp.queryAdd("lapse_advec_coeff", lapse_advec_coeff);
            if (std::min(std::abs(lapse_advec_coeff),
                         std::abs(lapse_advec_coeff - amrex::Real(1.0))) >
                std::numeric_limits<amrex::Real>::epsilon())
            {
                gauge_pp.warning("lapse_advec_coeff",
                                 "usually set to 0.0 or 1.0");
            }

            amrex::Real lapse_power = 1.;
            gauge_pp.queryAdd("lapse_power", lapse_power);
            if (std::abs(lapse_power - 1.0) >
                std::numeric_limits<amrex::Real>::epsilon())
            {
                gauge_pp.warning("lapse_power", "set to 1.0 for 1+log slicing");
            }

            amrex::Real lapse_coeff = 2.;
            gauge_pp.queryAdd("lapse_coeff", lapse_coeff);
            if (std::abs(lapse_coeff - 2.0) >
                std::numeric_limits<amrex::Real>::epsilon())
            {
                gauge_pp.warning("lapse_coeff", "set to 2.0 for 1+log slicing");
            }

            // Shift Evolution
            amrex::Real shift_Gamma_coeff = 0.75;
            gauge_pp.queryAdd("shift_Gamma_coeff", shift_Gamma_coeff);
            if (std::abs(shift_Gamma_coeff - 0.75) >
                std::numeric_limits<amrex::Real>::epsilon())
            {
                gauge_pp.warning("shift_Gamma_coeff", "usually set to 0.75");
            }

            amrex::Real shift_advec_coeff = 0.0;
            gauge_pp.queryAdd("shift_advec_coeff", shift_advec_coeff);
            if (std::min(std::abs(shift_advec_coeff),
                         std::abs(shift_advec_coeff - amrex::Real(1.0))) >
                std::numeric_limits<amrex::Real>::epsilon())
            {
                gauge_pp.warning("shift_advec_coeff",
                                 "usually set to 0.0 or 1.0");
            }

            amrex::Real eta = 1.0;
            gauge_pp.queryAdd("eta", eta);
            if (eta < 0.1 || eta > 10)
            {
                gauge_pp.warning(
                    "eta",
                    "usually O(1/M_ADM) so typically O(1) in code units");
            }
        }

        void fill_params()
        {
            GRParmParse gauge_pp("gauge");
            // lapse params:
            gauge_pp.get("lapse_advec_coeff", lapse_advec_coeff);
            gauge_pp.get("lapse_power", lapse_power);
            gauge_pp.get("lapse_coeff", lapse_coeff);
            // shift params:
            gauge_pp.get("shift_Gamma_coeff", shift_Gamma_coeff);
            gauge_pp.get("shift_advec_coeff", shift_advec_coeff);
            gauge_pp.get("eta", eta);
        }
    };

  protected:
    params_t m_params{};
    deriv_t m_deriv;

  public:
    MovingPunctureGauge(double a_dx) : m_deriv(a_dx) { m_params.fill_params(); }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    // NOLINTBEGIN(bugprone-easily-swappable-parameters,
    // readability-convert-member-functions-to-static)
    rhs_gauge(const amrex::CellData<amrex::Real> &rhs_cell_data,
              const CCZ4Vars &vars, const amrex::Real &advec_lapse,
              const Tensor::Rank1 &advec_shift, const Tensor::Rank1 &advec_B,
              const Tensor::Rank1 &advec_Gamma) const
    // NOLINTEND(bugprone-easily-swappable-parameters,
    // readability-convert-member-functions-to-static)
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

    /// Calculate the gauge RHS using the fully accumulated Gamma RHS.
    /** All vacuum and matter contributions to rhs[c_Gamma1 + i] must be added
     * before this function is called so that the B-field RHS uses the complete
     * time derivative of the conformal connection functions.
     */
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    calculate_gauge_rhs(int ix, int iy, int iz,
                        const amrex::Array4<amrex::Real> &rhs,
                        const amrex::Array4<const amrex::Real> &state) const
    {
        const amrex::CellData<amrex::Real> &rhs_cell_data =
            rhs.cellData(ix, iy, iz);
        const amrex::CellData<const amrex::Real> &state_cell_data =
            state.cellData(ix, iy, iz);

        CCZ4Vars vars(state_cell_data);

        Tensor::Rank1 shift_vector(
            {vars.shift(0), vars.shift(1), vars.shift(2)});

        amrex::Real advec_lapse =
            m_deriv.advec_scalar(ix, iy, iz, state, shift_vector, c_lapse);

        Tensor::Rank1 advec_shift =
            m_deriv.advec_vector(ix, iy, iz, state, shift_vector, c_shift1);

        Tensor::Rank1 advec_B =
            m_deriv.advec_vector(ix, iy, iz, state, shift_vector, c_B1);

        Tensor::Rank1 advec_Gamma =
            m_deriv.advec_vector(ix, iy, iz, state, shift_vector, c_Gamma1);

        rhs_gauge(rhs_cell_data, vars, advec_lapse, advec_shift, advec_B,
                  advec_Gamma);
    }
};

#endif /* MOVINGPUNCTUREGAUGE_HPP_ */
