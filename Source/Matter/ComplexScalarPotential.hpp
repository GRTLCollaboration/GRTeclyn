#ifndef COMPLEXSCALAR_POTENTIAL_HPP_
#define COMPLEXSCALAR_POTENTIAL_HPP_

#include <AMReX_REAL.H>
#include <cmath>

//! Shared complex-scalar potential V = 1/2 (m |Phi|)^2 - (lambda/4) |Phi|^4
//! on |Phi|^2 = phi1^2 + phi2^2.
class ComplexScalarPotential
{
  public:
    explicit ComplexScalarPotential(double a_mass = 0.0, double a_lambda = 0.0)
        : m_mass(a_mass), m_lambda(a_lambda)
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

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void compute_potential(
        amrex::Real &V_of_phi, amrex::Real &dVdphi1, amrex::Real &dVdphi2,
        amrex::Real phi1, amrex::Real phi2) const
    {
        const amrex::Real mod2 = phi1 * phi1 + phi2 * phi2;
        const amrex::Real mphi = m_mass * std::sqrt(mod2);
        const amrex::Real mod4 = mod2 * mod2;
        V_of_phi = 0.5 * mphi * mphi - 0.25 * m_lambda * mod4;
        if (mod2 <= 0.0)
        {
            dVdphi1 = 0.0;
            dVdphi2 = 0.0;
            return;
        }
        const amrex::Real dv_dmod2 = 0.5 * m_mass * m_mass - 0.5 * m_lambda * mod2;
        dVdphi1 = 2.0 * phi1 * dv_dmod2;
        dVdphi2 = 2.0 * phi2 * dv_dmod2;
    }

  private:
    double m_mass{};
    double m_lambda{};
};

#endif /* COMPLEXSCALAR_POTENTIAL_HPP_ */
