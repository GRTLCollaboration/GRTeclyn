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
        amrex::ParmParse pp("klein_gordon");
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
        // cos^2(beta*t) = 0.5*(1+cos(2*beta*t))

        amrex::Real numerator =
            m_alpha * std::sin(m_alpha * (t + m_t0)) * std::cosh(beta * x);
        amrex::Real denominator =
            0.5 * m_alpha * m_alpha * (1.0 + std::cosh(2.0 * beta * x)) +
            0.5 * beta * beta * (1.0 + std::cos(2.0 * beta * t));

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

    // NOLINTBEGIN(bugprone-easily-swappable-parameters
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    calculate_time_derivative(const amrex::Real x, const amrex::Real y,
                              const amrex::Real z, const amrex::Real t) const
    // NOLINTEND(bugprone-easily-swappable-parameters
    {
        // Sine Gordon 3D psuedo-breather solution

        amrex::Real beta = std::sqrt(1.0 - m_alpha * m_alpha);

        // Define these repeated calculations for convenience
        amrex::Real beta_time  = beta * (t + m_t0);
        amrex::Real sine_time  = m_alpha * std::sin(beta_time) / beta;
        amrex::Real cos_time   = m_alpha * beta * beta * std::cos(beta_time);
        amrex::Real cos_2_time = 1.0 - std::cos(2.0 * beta_time);

        amrex::Real x_term = std::atan(sine_time / std::cosh(m_alpha * x));
        amrex::Real y_term = std::atan(sine_time / std::cosh(m_alpha * y));
        amrex::Real z_term = std::atan(sine_time / std::cosh(m_alpha * z));

        amrex::Real dxdt_term = cos_time * std::cosh(m_alpha * x);
        dxdt_term /= (0.5 * beta * beta * (1.0 + std::cosh(2.0 * m_alpha * x)) +
                      0.5 * m_alpha * m_alpha * cos_2_time);

        amrex::Real dydt_term = cos_time * std::cosh(m_alpha * y);
        dydt_term /= (0.5 * beta * beta * (1.0 + std::cosh(2.0 * m_alpha * y)) +
                      0.5 * m_alpha * m_alpha * cos_2_time);

        amrex::Real dzdt_term = cos_time * std::cosh(m_alpha * z);
        dzdt_term /= (0.5 * beta * beta * (1.0 + std::cosh(2.0 * m_alpha * z)) +
                      0.5 * m_alpha * m_alpha * cos_2_time);

        return 64.0 *
               (dxdt_term * y_term * z_term + x_term * dydt_term * z_term +
                x_term * y_term * dzdt_term);
    };

    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void static compute_potential(
        amrex::Real &V_of_phi, amrex::Real &dVdphi, const amrex::Real &phi)
    // NOLINTEND(bugprone-easily-swappable-parameters)
    {
        V_of_phi = std::sin(phi);

        dVdphi = std::cos(phi);
    }
};
#endif // SINEGORDON_HPP_
