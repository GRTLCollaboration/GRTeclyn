#ifndef WAVE_HPP_
#define WAVE_HPP_

// C++ std lib includes
#include <cmath>
// AMReX includes
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
// KleinGordon includes
#include "StateVariables.hpp"

class Wave
{
  public:
    amrex::Real m_k_r{1.0};
    amrex::Real m_mass{0.0};
    amrex::Real m_t0{0.0};

    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
    Wave()
    {
        amrex::ParmParse pp("klein_gordon");
        pp.query("wave_vector", m_k_r);
        pp.query("scalar_mass", m_mass);
        pp.query("initial_time", m_t0);
    };

    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    calculate(const amrex::Real x, const amrex::Real y, const amrex::Real z,
              const amrex::Real t) const
    // NOLINTEND(bugprone-easily-swappable-parameters)
    {
        amrex::Real omega = m_k_r;

        // for the wave to be at the center of the grid, need to pass in
        // (x-x_midpt), (y-y_midpt) and (z-z_midpt)
        amrex::Real rr2 = x * x + y * y + z * z; // this is the radius

        return std::cos(m_k_r * rr2 - omega * (t + m_t0));
    };

    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    calculate_time_derivative(const amrex::Real x, const amrex::Real y,
                              const amrex::Real z, const amrex::Real t) const
    // NOLINTEND(bugprone-easily-swappable-parameters)
    {
        amrex::Real omega = m_k_r;

        // for the wave to be at the center of the grid, need to pass in
        // (x-x_midpt), (y-y_midpt) and (z-z_midpt)
        amrex::Real rr2 = x * x + y * y + z * z; // this is the radius

        return omega * std::sin(m_k_r * rr2 - omega * (t + m_t0));
    };

    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    template <class data_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    compute_potential(data_t &V_of_phi, data_t &dVdphi,
                      const amrex::Real &phi) const
    // NOLINTEND(bugprone-easily-swappable-parameters)
    {

        V_of_phi = 0.5 * m_mass * m_mass * phi * phi;

        dVdphi = m_mass * m_mass * phi;
    }
};

#endif // WAVE_HPP_
