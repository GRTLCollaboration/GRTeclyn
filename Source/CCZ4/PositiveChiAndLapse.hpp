/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// This compute class enforces the positive chi and lapse condition
#ifndef POSITIVECHIANDLAPSE_HPP_
#define POSITIVECHIANDLAPSE_HPP_

#include "Cell.hpp"
#include "StateVariables.hpp"

class PositiveChiAndLapse
{
  private:
    double m_min_chi;
    double m_min_lapse;

  public:
    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    //! Constructor for class
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE PositiveChiAndLapse(
        const double a_min_chi = 1e-4, const double a_min_lapse = 1e-4)
        : m_min_chi(a_min_chi), m_min_lapse(a_min_lapse)
    {
    }
    // NOLINTEND(bugprone-easily-swappable-parameters)

    template <class data_t>
    AMREX_GPU_HOST_DEVICE void
    operator()(const amrex::CellData<data_t> &cell) const
    {
        auto chi   = cell[c_chi];
        auto lapse = cell[c_lapse];

        chi   = (chi < m_min_chi) ? m_min_chi : chi;
        lapse = (lapse < m_min_lapse) ? m_min_lapse : chi;

        cell[c_chi]   = chi;
        cell[c_lapse] = lapse;
    }
};

#endif /* POSITIVECHIANDLAPSE_HPP_ */
