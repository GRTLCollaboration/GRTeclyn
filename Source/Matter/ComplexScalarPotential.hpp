#ifndef COMPLEXSCALAR_POTENTIAL_HPP_
#define COMPLEXSCALAR_POTENTIAL_HPP_

#include <AMReX_REAL.H>
#include <cmath>

//! Shared complex-scalar potential
//!   V = 1/2 (m |Phi|)^2 - (lambda/4) |Phi|^4 + (mu/6) |Phi|^6
//! on |Phi|^2 = phi1^2 + phi2^2.  The attractive quartic (lambda > 0) provides
//! self-binding; the sextic stabiliser (mu > 0) is REQUIRED for a stable 3D
//! Q-ball -- a pure attractive quartic is the critical case and either
//! collapses or disperses.  lambda = mu = 0 recovers the free massive field.
class ComplexScalarPotential
{
  public:
    explicit ComplexScalarPotential(double a_mass = 0.0, double a_lambda = 0.0,
                                    double a_mu = 0.0)
        : m_mass(a_mass), m_lambda(a_lambda), m_mu(a_mu)
    {
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE double mass() const
    {
        return m_mass;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE double lambda() const
    {
        return m_lambda;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE double mu() const
    {
        return m_mu;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void compute_potential(
        amrex::Real &V_of_phi, amrex::Real &dVdphi1, amrex::Real &dVdphi2,
        amrex::Real phi1, amrex::Real phi2) const
    {
        const amrex::Real mod2 = phi1 * phi1 + phi2 * phi2;
        const amrex::Real mod4 = mod2 * mod2;
        const amrex::Real mod6 = mod4 * mod2;
        // V = 1/2 m^2 |Phi|^2 - 1/4 lambda |Phi|^4 + 1/6 mu |Phi|^6
        V_of_phi = 0.5 * m_mass * m_mass * mod2 - 0.25 * m_lambda * mod4 +
                   (1.0 / 6.0) * m_mu * mod6;
        if (mod2 <= 0.0)
        {
            dVdphi1 = 0.0;
            dVdphi2 = 0.0;
            return;
        }
        // dV/d|Phi|^2 = 1/2 m^2 - 1/2 lambda |Phi|^2 + 1/2 mu |Phi|^4
        const amrex::Real dv_dmod2 = 0.5 * m_mass * m_mass -
                                     0.5 * m_lambda * mod2 +
                                     0.5 * m_mu * mod4;
        dVdphi1 = 2.0 * phi1 * dv_dmod2;
        dVdphi2 = 2.0 * phi2 * dv_dmod2;
    }

  private:
    double m_mass{};
    double m_lambda{};
    double m_mu{};
};

#endif /* COMPLEXSCALAR_POTENTIAL_HPP_ */
