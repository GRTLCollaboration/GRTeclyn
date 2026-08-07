/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef INTEGRATEDMOVINGPUNCTUREGAUGE_HPP_
#define INTEGRATEDMOVINGPUNCTUREGAUGE_HPP_

#include "CCZ4Vars.hpp"
#include "Cell.hpp"
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
class IntegratedMovingPunctureGauge
{
  public:
    using params_t = MovingPunctureGauge::params_t;

  protected:
    params_t m_params;

    /// Vars needed internally in 'compute'
    template <class data_t> struct Vars
    {
        Tensor::Rank1 shift;
        Tensor::Rank1 Gamma; //!< Conformal connection functions

        /// Defines the mapping between members of Vars and Chombo grid
        /// variables (enum in User_Variables)
        template <typename mapping_function_t>
        AMREX_GPU_HOST_DEVICE void
        enum_mapping(mapping_function_t mapping_function)
        {
            VarsTools::define_enum_mapping(
                mapping_function, GRInterval<c_shift1, c_shift3>(), shift);
            VarsTools::define_enum_mapping(
                mapping_function, GRInterval<c_Gamma1, c_Gamma3>(), Gamma);
        }
    };

  public:
    IntegratedMovingPunctureGauge() { m_params.fill_params(); }

    // set the initial B^i to the initial condition equivalent to:
    // \partial_t shift - advec_coeff * advec.shift = 0
    // Include in your Example in GRAMRLevel::initial_data as:
    // fillAllGhosts();
    // BoxLoops::loop(IntegratedMovingPunctureGauge(m_p.ccz4_params),
    // m_state_new, m_state_new, EXCLUDE_GHOST_CELLS);
    void compute(Cell<amrex::Real> current_cell) const
    {
        // TODO: Port this class
        // We've just removed templating over data_t
        const auto vars = current_cell.template load_vars<Vars>();

        Tensor::Rank1 B; // NOLINT(readability-identifier-length)
        FOR (i)
        {
            B(i) = m_params.shift_Gamma_coeff * vars.Gamma(i) -
                   m_params.eta * vars.shift(i);
        }

        current_cell.store_vars(B, GRInterval<c_B1, c_B3>());
    }

    template <template <class> class vars_t, class d2_vars_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    rhs_gauge(vars_t<amrex::Real> &rhs, const vars_t &vars<amrex::Real>,
              const vars_t<Tensor::Rank1> &d1, const d2_vars_t &d2,
              const vars_t<amrex::Real> &advec) const
    {
        rhs.lapse = m_params.lapse_advec_coeff * advec.lapse -
                    m_params.lapse_coeff *
                        pow(vars.lapse, m_params.lapse_power) *
                        (vars.K - 2 * vars.Theta);
        FOR (i)
        {
            rhs.shift(i) = m_params.shift_advec_coeff * advec.shift(i) +
                           m_params.shift_Gamma_coeff * vars.Gamma(i) -
                           m_params.eta * vars.shift(i) - vars.B(i);
            rhs.B(i) = 0.; // static, stays the same to save initial condition
        }
    }
};

#endif /* INTEGRATEDMOVINGPUNCTUREGAUGE_HPP_ */
