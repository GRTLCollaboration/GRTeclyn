/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// This compute class enforces the positive chi and alpha condition
#ifndef POSITIVECHIANDALPHA_HPP_
#define POSITIVECHIANDALPHA_HPP_

#include "Cell.hpp"
#include "StateVariables.hpp"
#include "simd.hpp"

class PositiveChiAndAlpha
{
  private:
    double m_min_chi;
    double m_min_lapse;

  public:
    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    //! Constructor for class
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE PositiveChiAndAlpha(
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

        chi   = simd_max(chi, m_min_chi);
        lapse = simd_max(lapse, m_min_lapse);

        cell[c_chi]   = chi;
        cell[c_lapse] = lapse;
    }
};

#endif /* POSITIVECHIANDALPHA_HPP_ */
