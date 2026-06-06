/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef OSCILLONPOTENTIAL_HPP_
#define OSCILLONPOTENTIAL_HPP_

#include "ScalarFieldVars.hpp"
#include <AMReX_REAL.H>

/// Real-scalar oscillon potential (metastable, not a true complex Q-ball).
/// V(phi) = 1/2 m^2 phi^2 - lambda/4 phi^4 + mu/6 phi^6
class OscillonPotential
{
  protected:
    double m_mass_sq;
    double m_lambda;
    double m_mu;

  public:
    OscillonPotential(double mass = 0.1, double lambda = 0.05, double mu = 0.01)
        : m_mass_sq(mass * mass), m_lambda(lambda), m_mu(mu)
    {
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    compute_potential(amrex::Real &V_of_phi, amrex::Real &dVdphi,
                      const ScalarFieldVars &vars) const
    {
        const amrex::Real phi = vars.phi();
        const amrex::Real phi2 = phi * phi;
        const amrex::Real phi4 = phi2 * phi2;
        V_of_phi = 0.5 * m_mass_sq * phi2 - 0.25 * m_lambda * phi4 +
                   (1.0 / 6.0) * m_mu * phi4 * phi2;
        dVdphi = m_mass_sq * phi - m_lambda * phi2 * phi + m_mu * phi4 * phi;
    }
};

#endif /* OSCILLONPOTENTIAL_HPP_ */
