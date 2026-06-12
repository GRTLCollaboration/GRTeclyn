#ifndef GRTRESNA_SCALAR_POTENTIAL_HPP_
#define GRTRESNA_SCALAR_POTENTIAL_HPP_

#include <AMReX_REAL.H>

//! Shared-sum scalar potential V = 1/2 (m s)^2 - (lambda/4) s^4 on s = sum phi_k.
class GRTresnaScalarPotential
{
  public:
    explicit GRTresnaScalarPotential(double a_mass = 0.0, double a_lambda = 0.0)
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

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    compute_potential(amrex::Real &V_of_phi, amrex::Real &dVdphi,
                      amrex::Real phi_sum) const
    {
        const amrex::Real mphi = m_mass * phi_sum;
        const amrex::Real s2   = phi_sum * phi_sum;
        const amrex::Real s4   = s2 * s2;
        V_of_phi               = 0.5 * mphi * mphi - 0.25 * m_lambda * s4;
        dVdphi                 = m_mass * mphi - m_lambda * s2 * phi_sum;
    }

  private:
    double m_mass{};
    double m_lambda{};
};

#endif /* GRTRESNA_SCALAR_POTENTIAL_HPP_ */
