#ifndef GRTRESNA_SCALAR_POTENTIAL_HPP_
#define GRTRESNA_SCALAR_POTENTIAL_HPP_

#include <AMReX_REAL.H>

//! Massive scalar potential V = 1/2 m^2 phi^2 matching GRTresna scalar_mass.
class GRTresnaScalarPotential
{
  public:
    explicit GRTresnaScalarPotential(double a_mass = 0.0) : m_mass(a_mass) {}

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE double mass() const
    {
        return m_mass;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    compute_potential(amrex::Real &V_of_phi, amrex::Real &dVdphi,
                      amrex::Real phi_sum) const
    {
        const amrex::Real mphi = m_mass * phi_sum;
        V_of_phi               = 0.5 * mphi * mphi;
        dVdphi                 = m_mass * mphi;
    }

  private:
    double m_mass{};
};

#endif /* GRTRESNA_SCALAR_POTENTIAL_HPP_ */
