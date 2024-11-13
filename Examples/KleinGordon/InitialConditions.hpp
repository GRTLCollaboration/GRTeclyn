#ifndef INITIAL_CONDITIONS_H_
#define INITIAL_CONDITIONS_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <cmath>

class InitialConditions
{

  private:

    // this is either the initial value of alpha for Sine Gordon ICs or could
    // also be the wavenumber for Wave ICs

    amrex::Real m_initial{1};

  public:
    InitialConditions(const amrex::Real a_initial) : m_initial(a_initial){};

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
    breather_solution(const amrex::Real x, const amrex::Real t) const
    {
        // Sine Gordon 1D breather solution
        amrex::Real alpha = m_initial;
        amrex::Real beta  = std::sqrt(1.0 - alpha * alpha);

        return 4 * std::atan(beta * std::cos(alpha * t) / alpha /
                             std::cosh(beta * x));
    }

    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    breather_solution(const amrex::Real x, const amrex::Real y,
                      const amrex::Real z, const amrex::Real t) const
    // NOLINTEND(bugprone-easily-swappable-parameters)
    {
        amrex::Real alpha = m_initial;
        amrex::Real beta  = std::sqrt(1.0 - alpha * alpha);

        // First derivative of Sine Gordon 1D breather solution
        amrex::Real x1_origin = 0; // for two different breathers
        amrex::Real x2_origin = 0; // located at x1_origin and x2_origin
        amrex::Real v         = 0;

        amrex::Real y1_origin = t - v * x + x1_origin;
        amrex::Real y2_origin = x - v * t + x2_origin;

        amrex::Real numerator =
            alpha * std::sin(alpha * y1_origin) * std::cosh(beta * y2_origin);
        amrex::Real denominator = alpha * alpha * std::cosh(beta * y2_origin) *
                                      std::cosh(beta * y2_origin) +
                                  beta * beta * std::cos(alpha * y1_origin) *
                                      std::cos(alpha * y1_origin);

        return -4 * alpha * beta * numerator / denominator;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
    breather_solution_deriv(const amrex::Real x, const amrex::Real t) const
    {
        amrex::Real alpha = m_initial;
        amrex::Real beta  = std::sqrt(1.0 - alpha * alpha);

        // First derivative of Sine Gordon 1D breather solution
        amrex::Real x1_origin = 0;
        amrex::Real x2_origin = 0;
        amrex::Real v         = 0;

        amrex::Real y1_origin = t - v * x + x1_origin;
        amrex::Real y2_origin = x - v * t + x2_origin;

        amrex::Real numerator =
            alpha * std::sin(alpha * y1_origin) * std::cosh(beta * y2_origin);
        amrex::Real denominator = alpha * alpha * std::cosh(beta * y2_origin) *
                                      std::cosh(beta * y2_origin) +
                                  beta * beta * std::cos(alpha * y1_origin) *
                                      std::cos(alpha * y1_origin);

        return -4 * alpha * beta * numerator / denominator;
    }

    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    travelling_wave(const amrex::Real x, const amrex::Real y,
                    const amrex::Real z, const amrex::Real t) const
    // NOLINTEND(bugprone-easily-swappable-parameters)
    {
        amrex::Real k_r   = m_initial;
        amrex::Real omega = k_r;

        // for the wave to be at the center of the grid, need to pass in
        // (x-x_midpt), (y-y_midpt) and (z-z_midpt)
        amrex::Real rr2 = x * x + y * y + z * z; // this is the radius

        return std::cos(k_r * rr2 - omega * t);
    }

    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    travelling_wave_deriv(const amrex::Real x, const amrex::Real y,
                          const amrex::Real z, const amrex::Real t) const
    // NOLINTEND(bugprone-easily-swappable-parameters)
    {
        amrex::Real k_r   = m_initial;
        amrex::Real omega = k_r;

        // for the wave to be at the center of the grid, need to pass in
        // (x-x_midpt), (y-y_midpt) and (z-z_midpt)
        amrex::Real rr2 = x * x + y * y + z * z; // this is the radius

        return omega * std::sin(k_r * rr2 - omega * t);
    }
};

#endif /* INITIAL_CONDITONS_H_ */
