/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SCALARFIELDVARS_HPP_
#define SCALARFIELDVARS_HPP_

#include "CCZ4Vars.hpp"
#include "StateVariables.hpp"
#include "Tensor.hpp"

class ScalarFieldVars : public CCZ4Vars
{
  public:
    AMREX_GPU_DEVICE inline ScalarFieldVars(
        const amrex::CellData<amrex::Real> &input_cell_data)
        : CCZ4Vars(input_cell_data)
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

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void store_phi(amrex::Real phi)
    {
        cell_data[c_phi] = phi;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void store_Pi(amrex::Real Pi)
    {
        cell_data[c_Pi] = Pi;
    }
};

#endif /* SCALARFIELDVARS_HPP */