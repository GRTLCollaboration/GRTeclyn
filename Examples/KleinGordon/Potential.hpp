#include "StateVariables.hpp"
#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>

#include <cmath>

class Potential
{
  private:
    amrex::Real m_mass;

  public:
    Potential(const amrex::Real mass) : m_mass(mass){};
    ~Potential(){};

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    compute_phi_sq(int i, int j, int k,
                   const amrex::Array4<amrex::Real> &phi) const
    {
        auto phi2 = phi(i, j, k, c_phi) * phi(i, j, k, c_phi);

        phi(i, j, k, c_phi) += 0.5 * m_mass * m_mass * phi2;
    }
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    compute_sine_gordon(int i, int j, int k,
                        const amrex::Array4<amrex::Real> &phi) const
    {

        auto phi2 = phi(i, j, k, c_phi) * phi(i, j, k, c_phi);

        phi(i, j, k, c_phi) += std::sin(phi2);
    }
};
