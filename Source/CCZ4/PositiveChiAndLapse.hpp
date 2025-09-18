/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// This compute class enforces the positive chi and lapse condition
#ifndef POSITIVECHIANDLAPSE_HPP_
#define POSITIVECHIANDLAPSE_HPP_

#include "CCZ4Vars.hpp"
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

    AMREX_GPU_DEVICE void
    operator()(int ix, int iy, int iz,
               const amrex::Array4<amrex::Real> &state) const
    {

        const amrex::CellData<amrex::Real> &state_cell_data =
            state.cellData(ix, iy, iz);
        CCZ4Vars vars(state_cell_data);

        amrex::Real chi   = vars.chi();
        amrex::Real lapse = vars.lapse();

        state_cell_data[c_chi]   = std::max(chi, m_min_chi);
        state_cell_data[c_lapse] = std::max(lapse, m_min_lapse);
    }
};

#endif /* POSITIVECHIANDLAPSE_HPP_ */
