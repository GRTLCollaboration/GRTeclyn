/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef COMPLEXSCALARFIELDVARS_HPP_
#define COMPLEXSCALARFIELDVARS_HPP_

#include "CCZ4Vars.hpp"
#include "StateVariables.hpp"
#include "Tensor.hpp"

//! State accessors for a complex scalar stored as two real components
//! (phi1, Pi1) and (phi2, Pi2). The complex field is Phi = phi1 + i phi2.
class ComplexScalarFieldVars : public CCZ4Vars
{
  public:
    AMREX_GPU_DEVICE
    ComplexScalarFieldVars(
        const amrex::CellData<const amrex::Real> &input_cell_data)
        : CCZ4Vars(input_cell_data)
    {
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &phi1() const
    {
        return cell_data[c_phi];
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &Pi1() const
    {
        return cell_data[c_Pi];
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &phi2() const
    {
        return cell_data[c_phi2];
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &Pi2() const
    {
        return cell_data[c_Pi2];
    }
};

#endif /* COMPLEXSCALARFIELDVARS_HPP_ */
