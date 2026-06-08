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
        compute_potential_value(V_of_phi, dVdphi, vars.phi());
    }

    //! Value-based overload used by the complex scalar field, which evaluates
    //! the potential separately on each real component. The quadratic phantom
    //! potential is separable: V(|Phi|^2) = sum_a 0.5 m^2 phi_a^2, so the
    //! per-component evaluation is exact.
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
    compute_potential_value(amrex::Real &V_of_phi, amrex::Real &dVdphi,
                            amrex::Real phi) const
    {
        V_of_phi = 0.5 * m_mass_sq * phi * phi;
        dVdphi = m_mass_sq * phi;
    }
};

#endif /* PHANTOMDECAYPOTENTIAL_HPP_ */
