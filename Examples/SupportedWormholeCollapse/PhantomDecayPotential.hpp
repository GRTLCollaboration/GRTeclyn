/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef PHANTOMDECAYPOTENTIAL_HPP_
#define PHANTOMDECAYPOTENTIAL_HPP_

#include "ScalarFieldVars.hpp"
#include "simd.hpp"
#include "Tensor.hpp"

class PhantomDecayPotential
{
  protected:
    double m_mass_sq;

  public:
    PhantomDecayPotential(double mass = 0.0) : m_mass_sq(mass * mass) {}

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
    compute_potential(amrex::Real &V_of_phi, amrex::Real &dVdphi,
                      const ScalarFieldVars &vars) const
    {
        const amrex::Real phi = vars.phi();

        V_of_phi = 0.5 * m_mass_sq * phi * phi;
        dVdphi = m_mass_sq * phi;
    }
};

#endif /* PHANTOMDECAYPOTENTIAL_HPP_ */
