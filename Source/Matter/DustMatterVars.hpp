/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef DUSTMATTERVARS_HPP_
#define DUSTMATTERVARS_HPP_

#include "CCZ4Vars.hpp"

class DustMatterVars : public CCZ4Vars
{
  public:
    AMREX_GPU_DEVICE
    DustMatterVars(const amrex::CellData<const amrex::Real> &input_cell_data)
        : CCZ4Vars(input_cell_data)
    {
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &
    dust_rho() const
    {
        return cell_data[c_dust_rho];
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &
    dust_v(int i) const
    {
        return cell_data[c_dust_v1 + i];
    }
};

#endif /* DUSTMATTERVARS_HPP_ */
