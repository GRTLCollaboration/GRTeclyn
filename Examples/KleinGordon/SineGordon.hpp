/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SINEGORDON_HPP_
#define SINEGORDON_HPP_

// C++ std lib includes
#include <cmath>
// AMReX includes
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
// KleinGordon includes
#include "StateVariables.hpp"
#include "GRParmParse.hpp"

class SineGordon
{
  public:

    struct params_t
    {
        amrex::Real alpha{};
        static void check_params()
        {
            // These are parameters specfic to the Sine Gordon example
            GRParmParse sg_pp("sine_gordon");
            amrex::Real alpha{1.0};
            sg_pp.queryAdd("alpha", alpha);
        }

        void fill_params()
        {
            GRParmParse sg_pp("sine_gordon");
            sg_pp.get("alpha", alpha);
        }
    };

    amrex::Real m_t0{};
    params_t m_params;

    SineGordon()
    {
        GRParmParse pp;
        pp.query("klein_gordon.initial_time", m_t0);

        m_params.fill_params();
    };

    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    calculate(const amrex::Real x, const amrex::Real t) const
    // NOLINTEND(bugprone-easily-swappable-parameters)
    {
        // Sine Gordon 1D breather solution
        amrex::Real beta = std::sqrt(1.0 - m_params.alpha * m_params.alpha);

        return 4.0 * std::atan(beta * std::cos(m_params.alpha * (t + m_t0)) / m_params.alpha /
                               std::cosh(beta * x));
    };
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
    calculate_time_derivative(const amrex::Real x, const amrex::Real t) const
    {
        amrex::Real beta = std::sqrt(1.0 - m_params.alpha * m_params.alpha);

        // First derivative of Sine Gordon 1D breather solution

        // Using cosh^2 (beta*x) = 1/2(1+cosh(2*beta*x)) and
        // cos^2(beta*t) = 0.5*(1+cos(2*beta*t))

        amrex::Real numerator =
            m_params.alpha * std::sin(m_params.alpha * (t + m_t0)) * std::cosh(beta * x);
        amrex::Real denominator =
            0.5 * m_params.alpha * m_params.alpha * (1.0 + std::cosh(2.0 * beta * x)) +
            0.5 * beta * beta * (1.0 + std::cos(2.0 * beta * t));

        return -4.0 * m_params.alpha * beta * numerator / denominator;
    };

    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    calculate(const amrex::Real x, const amrex::Real y, const amrex::Real z,
              const amrex::Real t) const
    // NOLINTEND(bugprone-easily-swappable-parameters)
    {
        amrex::Real beta = std::sqrt(1.0 - m_params.alpha * m_params.alpha);

        amrex::Real numerator = m_params.alpha * std::sin(beta * (t + m_t0)) / beta;

        // Sine Gordon 3D psuedo-breather solution

        return 4.0 * 4.0 * 4.0 * std::atan(numerator / std::cosh(m_params.alpha * x)) *
               std::atan(numerator / std::cosh(m_params.alpha * y)) *
               std::atan(numerator / std::cosh(m_params.alpha * z));
    };

    // NOLINTBEGIN(bugprone-easily-swappable-parameters
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    calculate_time_derivative(const amrex::Real x, const amrex::Real y,
                              const amrex::Real z, const amrex::Real t) const
    // NOLINTEND(bugprone-easily-swappable-parameters
    {
        // Sine Gordon 3D psuedo-breather solution

        amrex::Real beta = std::sqrt(1.0 - m_params.alpha * m_params.alpha);

        // Define these repeated calculations for convenience
        amrex::Real beta_time  = beta * (t + m_t0);
        amrex::Real sine_time  = m_params.alpha * std::sin(beta_time) / beta;
        amrex::Real cos_time   = m_params.alpha * beta * beta * std::cos(beta_time);
        amrex::Real cos_2_time = 1.0 - std::cos(2.0 * beta_time);

        amrex::Real x_term = std::atan(sine_time / std::cosh(m_params.alpha * x));
        amrex::Real y_term = std::atan(sine_time / std::cosh(m_params.alpha * y));
        amrex::Real z_term = std::atan(sine_time / std::cosh(m_params.alpha * z));

        amrex::Real dxdt_term = cos_time * std::cosh(m_params.alpha * x);
        dxdt_term /= (0.5 * beta * beta * (1.0 + std::cosh(2.0 * m_params.alpha * x)) +
                      0.5 * m_params.alpha * m_params.alpha * cos_2_time);

        amrex::Real dydt_term = cos_time * std::cosh(m_params.alpha * y);
        dydt_term /= (0.5 * beta * beta * (1.0 + std::cosh(2.0 * m_params.alpha * y)) +
                      0.5 * m_params.alpha * m_params.alpha * cos_2_time);

        amrex::Real dzdt_term = cos_time * std::cosh(m_params.alpha * z);
        dzdt_term /= (0.5 * beta * beta * (1.0 + std::cosh(2.0 * m_params.alpha * z)) +
                      0.5 * m_params.alpha * m_params.alpha * cos_2_time);

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
