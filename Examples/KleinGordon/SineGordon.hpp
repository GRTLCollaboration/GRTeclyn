#ifndef SINEGORDON_HPP_
#define SINEGORDON_HPP_

// C++ std lib includes
#include <cmath>
// AMReX includes
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
// KleinGordon includes
#include "StateVariables.hpp"

class SineGordon
{
  public:
    amrex::Real m_alpha{0.7};
    amrex::Real m_t0{0.0};

    SineGordon()
    {
        amrex::ParmParse pp;
        pp.query("alpha", m_alpha);
        pp.query("initial_time", m_t0);
    };

    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    calculate(const amrex::Real x, const amrex::Real t) const
    // NOLINTEND(bugprone-easily-swappable-parameters)
    {
        // Sine Gordon 1D breather solution
        amrex::Real beta = std::sqrt(1.0 - m_alpha * m_alpha);

        return 4.0 * std::atan(beta * std::cos(m_alpha * (t + m_t0)) / m_alpha /
                               std::cosh(beta * x));
    };
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
    calculate_time_derivative(const amrex::Real x, const amrex::Real t) const
    {
        amrex::Real beta = std::sqrt(1.0 - m_alpha * m_alpha);

        // First derivative of Sine Gordon 1D breather solution

        // Using cosh^2 (beta*x) = 1/2(1+cosh(2*beta*x)) and
        // cos^2(alpha*x) = 0.5*(1+cos(2*beta*x))

        amrex::Real numerator =
            m_alpha * std::sin(m_alpha * (t + m_t0)) * std::cosh(beta * x);
        amrex::Real denominator =
            0.5 * m_alpha * m_alpha * (1.0 + std::cosh(2.0 * beta * x)) +
            0.5 * beta * beta * (1.0 + std::cos(2.0 * beta * x));

        return -4.0 * m_alpha * beta * numerator / denominator;
    };

    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    calculate(const amrex::Real x, const amrex::Real y, const amrex::Real z,
              const amrex::Real t) const
    // NOLINTEND(bugprone-easily-swappable-parameters)
    {
        amrex::Real beta = std::sqrt(1.0 - m_alpha * m_alpha);

        amrex::Real numerator = m_alpha * std::sin(beta * (t + m_t0)) / beta;

        // Sine Gordon 3D psuedo-breather solution

        return 4.0 * 4.0 * 4.0 * std::atan(numerator / std::cosh(m_alpha * x)) *
               std::atan(numerator / std::cosh(m_alpha * y)) *
               std::atan(numerator / std::cosh(m_alpha * z));
    };

    // NOLINTBEGIN(bugprone-easily-swappable-parameters,
    // readability-convert-member-functions-to-static)
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    calculate_time_derivative(const amrex::Real x, const amrex::Real y,
                              const amrex::Real z, const amrex::Real t) const
    // NOLINTEND(bugprone-easily-swappable-parameters,
    // readability-convert-member-functions-to-static)
    {
        // Sine Gordon 3D psuedo-breather solution

        return 0.0;
    };

    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    template <class data_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    compute_potential(data_t &V_of_phi, data_t &dVdphi,
                      const amrex::Real &phi) const
    // NOLINTEND(bugprone-easily-swappable-parameters)
    {
        V_of_phi = std::sin(phi);

        dVdphi = std::cos(phi);
    }
};
#endif // SINEGORDON_HPP_
