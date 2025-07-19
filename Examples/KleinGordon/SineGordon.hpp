#ifndef SINEGORDON_HPP_
#define SINEGORDON_HPP_

// C++ std lib includes
#include <cmath>
// AMReX includes
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
// GRTeclyn includes
#include "simd.hpp"
// KleinGordon includes
#include "StateVariables.hpp"

class SineGordon
{
  public:
    amrex::Real m_alpha{0.7};
    amrex::Real m_t0{0.0};

    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
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

        return 4 * std::atan(beta * std::cos(m_alpha * (t + m_t0)) / m_alpha /
                             std::cosh(beta * x));
    };
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
    calculate_derivative(const amrex::Real x, const amrex::Real t) const
    {
        amrex::Real beta = std::sqrt(1.0 - m_alpha * m_alpha);

        // First derivative of Sine Gordon 1D breather solution
        amrex::Real x1_origin = 0;
        amrex::Real x2_origin = 0;
        amrex::Real v         = 0;

        amrex::Real y1_origin = (t + m_t0) - v * x + x1_origin;
        amrex::Real y2_origin = x - v * (t + m_t0) + x2_origin;

        amrex::Real numerator = m_alpha * std::sin(m_alpha * y1_origin) *
                                std::cosh(beta * y2_origin);
        amrex::Real denominator = m_alpha * m_alpha *
                                      std::cosh(beta * y2_origin) *
                                      std::cosh(beta * y2_origin) +
                                  beta * beta * std::cos(m_alpha * y1_origin) *
                                      std::cos(m_alpha * y1_origin);

        return -4 * m_alpha * beta * numerator / denominator;
    };

    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    calculate(const amrex::Real x, const amrex::Real y, const amrex::Real z,
              const amrex::Real t) const
    // NOLINTEND(bugprone-easily-swappable-parameters)
    {
        amrex::Real beta = std::sqrt(1.0 - m_alpha * m_alpha);

        // Sine Gordon 3D psuedo-breather solution

        return 4 * 4 * 4 *
               std::atan(m_alpha * std::sin(beta * (t + m_t0)) / beta /
                         std::cosh(m_alpha * x)) *
               std::atan(m_alpha * std::sin(beta * (t + m_t0)) / beta /
                         std::cosh(m_alpha * y)) *
               std::atan(m_alpha * std::sin(beta * (t + m_t0)) / beta /
                         std::cosh(m_alpha * z));
    };

    // NOLINTBEGIN(bugprone-easily-swappable-parameters,
    // readability-convert-member-functions-to-static)
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    calculate_derivative(const amrex::Real x, const amrex::Real y,
                         const amrex::Real z, const amrex::Real t) const
    // NOLINTEND(bugprone-easily-swappable-parameters,
    // readability-convert-member-functions-to-static)
    {
        // Sine Gordon 3D psuedo-breather solution

        return 0.0;
    };

    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    template <class data_t, template <typename> class vars_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    compute_potential(data_t &V_of_phi, data_t &dVdphi,
                      const vars_t<data_t> &vars) const
    // NOLINTEND(bugprone-easily-swappable-parameters)
    {
        V_of_phi = std::sin(vars.phi);

        dVdphi = std::cos(vars.phi);
    }
};
#endif // SINEGORDON_HPP_
