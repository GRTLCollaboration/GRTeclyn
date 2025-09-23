/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef CONSTSCALARFIELDVARS_HPP_
#define CONSTSCALARFIELDVARS_HPP_

#include "ConstCCZ4Vars.hpp"
#include "StateVariables.hpp"
#include "Tensor.hpp"

class ConstScalarFieldVars : public ConstCCZ4Vars
{
  public:
    AMREX_GPU_DEVICE ConstScalarFieldVars(
        const amrex::CellData<const amrex::Real> &input_cell_data)
        : ConstCCZ4Vars(input_cell_data)
    {
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &phi() const
    {
        return cell_data[c_phi];
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &Pi() const
    {
        return cell_data[c_Pi];
    }
};

#endif /* CONSTSCALARFIELDVARS_HPP */