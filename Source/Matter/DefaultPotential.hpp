/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef DEFAULTPOTENTIAL_HPP_
#define DEFAULTPOTENTIAL_HPP_

#include "ConstScalarFieldVars.hpp"
#include "Tensor.hpp"
#include <AMReX_REAL.H>

class DefaultPotential
{
  public:
    //! The constructor
    DefaultPotential() = default;

    //! Set the potential function for the scalar field here to zero
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
    compute_potential(amrex::Real &V_of_phi, amrex::Real &dVdphi,
                      const ConstScalarFieldVars &vars) const
    {
        // The potential value at phi
        V_of_phi = 0.0;

        // The potential gradient at phi
        dVdphi = 0.0;
    }
};

#endif /* DEFAULTPOTENTIAL_HPP_ */
