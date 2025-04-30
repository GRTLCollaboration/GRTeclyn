/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef MOVINGPUNCTUREGAUGEWITHMATTER_HPP_
#define MOVINGPUNCTUREGAUGEWITHMATTER_HPP_

#include "CCZ4RHSWithMatter.hpp"
#include "DimensionDefinitions.hpp"
// #include "EMTensor.hpp"
#include "Tensor.hpp"

#include <cmath>

/**
 * This class implements a slightly more generic version of the moving puncture
 * gauge. In particular it uses a Bona-Masso slicing condition of the form
 * f(lapse) = -c*lapse^(p-2)
 * and a Gamma-driver shift condition
 **/

/// This class adds the matter terms to the RHS of the gauge equation
/// for the moving puncture gauge

class MovingPunctureGaugeWithMatter : public MovingPunctureGauge
{

  public:
    MovingPunctureGaugeWithMatter(const params_t &a_params)
        : MovingPunctureGauge(a_params)
    {
    }

    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    template <class data_t, template <typename> class vars_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void rhs_gauge_add_matter_terms(
        vars_t<data_t> &matter_rhs, const vars_t<data_t> &matter_vars,
        Tensor<2, data_t, 3> h_UU, const emtensor_t<data_t> emtensor,
        const double G_Newton) const
    // NOLINTEND(bugprone-easily-swappable-parameters)
    {
        FOR (i)
        {
            data_t matter_term_Gamma = 0.0;
            FOR (j)
            {
                matter_term_Gamma += -16.0 * M_PI * G_Newton *
                                     matter_vars.lapse * h_UU[i][j] *
                                     emtensor.j[j];
            }

            matter_rhs.B[i] += matter_term_Gamma;
        }
    }
};

#endif /* MOVINGPUNCTUREGAUGEWITHMATTER_HPP_ */
