/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef MOVINGPUNCTUREGAUGEWITHMATTER_HPP_
#define MOVINGPUNCTUREGAUGEWITHMATTER_HPP_

#include "CCZ4Geometry.hpp"
#include "CCZ4Vars.hpp"
#include "DimensionDefinitions.hpp"
#include "MovingPunctureGauge.hpp"

#include <cmath>

/**
 * This class implements a slightly more generic version of the moving puncture
 * gauge. In particular it uses a Bona-Masso slicing condition of the form
 * f(lapse) = -c*lapse^(p-2)
 * and a Gamma-driver shift condition
 **/

/// This class adds the matter terms to the RHS of the gauge equation
/// for the moving puncture gauge

template <class deriv_t = FourthOrderDerivatives>
class MovingPunctureGaugeWithMatter : public MovingPunctureGauge<deriv_t>
{
    using base_t = MovingPunctureGauge<deriv_t>;

  public:
    MovingPunctureGaugeWithMatter(double a_dx) : base_t(a_dx) {}

    using base_t::calculate_gauge_rhs;

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    calculate_gauge_rhs(int ix, int iy, int iz,
                        const amrex::Array4<amrex::Real> &rhs,
                        const amrex::Array4<const amrex::Real> &state,
                        const Tensor::Rank2 &h_UU,
                        const einstein_sources_t &source) const
    {
        base_t::calculate_gauge_rhs(ix, iy, iz, rhs, state);

        const amrex::CellData<amrex::Real> &rhs_cell_data =
            rhs.cellData(ix, iy, iz);
        const amrex::CellData<const amrex::Real> &state_cell_data =
            state.cellData(ix, iy, iz);
        const CCZ4Vars vars(state_cell_data);

        rhs_gauge_add_matter_terms(rhs_cell_data, vars, h_UU, source);
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE static void
    rhs_gauge_add_matter_terms(const amrex::CellData<amrex::Real> &rhs,
                               const CCZ4Vars &vars, const Tensor::Rank2 &h_UU,
                               const einstein_sources_t &source)
    {
        FOR (i)
        {
            amrex::Real matter_term_Gamma = 0.0;
            FOR (j)
            {
                matter_term_Gamma -=
                    2.0 * vars.lapse() * h_UU(i, j) * source.j(j);
            }
            rhs[c_B1 + i] += matter_term_Gamma;
        }
    }
};

#endif /* MOVINGPUNCTUREGAUGEWITHMATTER_HPP_ */
