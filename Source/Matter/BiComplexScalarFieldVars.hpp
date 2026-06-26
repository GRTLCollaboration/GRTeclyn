/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef BICOMPLEXSCALARFIELDVARS_HPP_
#define BICOMPLEXSCALARFIELDVARS_HPP_

#include "CCZ4Vars.hpp"
#include "StateVariables.hpp"
#include "Tensor.hpp"

//! State accessors for TWO complex scalar fields stored as real components:
//! a CANONICAL field  Phi+ = phi1p + i phi2p  (conjugate momenta Pi1p, Pi2p)
//! and a PHANTOM field Phi- = phi1m + i phi2m (conjugate momenta Pi1m, Pi2m).
//! The phantom field couples to gravity with a flipped (negative) sign.
class BiComplexScalarFieldVars : public CCZ4Vars
{
  public:
    AMREX_GPU_DEVICE
    BiComplexScalarFieldVars(
        const amrex::CellData<const amrex::Real> &input_cell_data)
        : CCZ4Vars(input_cell_data)
    {
    }

    // Canonical (Phi+) field.
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &
    phi1p() const { return cell_data[c_phi]; }
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &
    Pi1p() const { return cell_data[c_Pi]; }
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &
    phi2p() const { return cell_data[c_phi2]; }
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &
    Pi2p() const { return cell_data[c_Pi2]; }

    // Phantom (Phi-) field.
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &
    phi1m() const { return cell_data[c_phi_m]; }
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &
    Pi1m() const { return cell_data[c_Pi_m]; }
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &
    phi2m() const { return cell_data[c_phi2_m]; }
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &
    Pi2m() const { return cell_data[c_Pi2_m]; }
};

#endif /* BICOMPLEXSCALARFIELDVARS_HPP_ */
