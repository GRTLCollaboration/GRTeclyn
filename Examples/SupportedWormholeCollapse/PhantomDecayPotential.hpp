/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef PHANTOMDECAYPOTENTIAL_HPP_
#define PHANTOMDECAYPOTENTIAL_HPP_

#include "ScalarFieldVars.hpp"
#include "simd.hpp"
#include "Tensor.hpp"

//! A simple mass-term potential for a scalar field: V(phi) = 1/2 m^2 phi^2
//! When used with ExoticScalarField (which flips the sign of the energy-momentum tensor),
//! this creates an effective tachyonic instability that causes the phantom field to naturally disperse.
class PhantomDecayPotential
{
  protected:
    double m_mass_sq;

  public:
    //! Constructor. Parameter 'mass' controls the decay timescale (~1/m).
    PhantomDecayPotential(double mass = 0.0) : m_mass_sq(mass * mass) {}

    //! Computes the potential V(phi) and its derivative dV/dphi
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
    compute_potential(amrex::Real &V_of_phi, amrex::Real &dVdphi,
                      const ScalarFieldVars &vars) const
    {
        const amrex::Real phi = vars.phi();

        // The potential value at phi: V(phi) = 1/2 m^2 phi^2
        V_of_phi = 0.5 * m_mass_sq * phi * phi;

        // The potential gradient at phi: dV/dphi = m^2 phi
        dVdphi = m_mass_sq * phi;
    }
};

#endif /* PHANTOMDECAYPOTENTIAL_HPP_ */
