#ifndef INITIAL_CONDITIONS_H_
#define INITIAL_CONDITIONS_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <cmath>

class InitialConditions
{

  private:
    const amrex::Real m_alpha;
    const amrex::Real m_beta;

    const amrex::Real m_k_r;
    const amrex::Real m_omega;

  public:
    InitialConditions(const amrex::Real alpha, const amrex::Real k_r)
        : m_alpha(alpha), m_beta(std::sqrt(1.0 - alpha * alpha)), m_k_r(k_r),
          m_omega(k_r){};

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    breather_solution(const amrex::Real x, const amrex::Real t) const
    {
        // Sine Gordon 1D breather solution

        return 4 * std::atan(m_beta * std::cos(m_alpha * t) / m_alpha /
                             std::cosh(m_beta * x));
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    breather_solution(const amrex::Real x, const amrex::Real y,
                      const amrex::Real z, const amrex::Real t) const
    {
        // First derivative of Sine Gordon 1D breather solution
        amrex::Real x1 = 0;
        amrex::Real x2 = 0;
        amrex::Real v  = 0;

        amrex::Real y1 = t - v * x + x1;
        amrex::Real y2 = x - v * t + x2;

        amrex::Real numerator =
            m_alpha * std::sin(m_alpha * y1) * std::cosh(m_beta * y2);
        amrex::Real denominator =
            m_alpha * m_alpha * std::cosh(m_beta * y2) *
                std::cosh(m_beta * y2) +
            m_beta * m_beta * std::cos(m_alpha * y1) * std::cos(m_alpha * y1);

        return -4 * m_alpha * m_beta * numerator / denominator;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    breather_solution_deriv(const amrex::Real x, const amrex::Real t) const
    {
        // First derivative of Sine Gordon 1D breather solution
        amrex::Real x1 = 0;
        amrex::Real x2 = 0;
        amrex::Real v  = 0;

        amrex::Real y1 = t - v * x + x1;
        amrex::Real y2 = x - v * t + x2;

        amrex::Real numerator =
            m_alpha * std::sin(m_alpha * y1) * std::cosh(m_beta * y2);
        amrex::Real denominator =
            m_alpha * m_alpha * std::cosh(m_beta * y2) *
                std::cosh(m_beta * y2) +
            m_beta * m_beta * std::cos(m_alpha * y1) * std::cos(m_alpha * y1);

        return -4 * m_alpha * m_beta * numerator / denominator;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    travelling_wave(const amrex::Real x, const amrex::Real y,
                    const amrex::Real z, const amrex::Real t) const
    {
        // for the wave to be at the center of the grid, need to pass in
        // (x-x_midpt), (y-y_midpt) and (z-z_midpt)
        amrex::Real rr2 = x * x + y * y + z * z; // this is the radius

        return std::cos(m_k_r * rr2 - m_omega * t);
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    travelling_wave_deriv(const amrex::Real x, const amrex::Real y,
                          const amrex::Real z, const amrex::Real t) const
    {
        // for the wave to be at the center of the grid, need to pass in
        // (x-x_midpt), (y-y_midpt) and (z-z_midpt)
        amrex::Real rr2 = x * x + y * y + z * z; // this is the radius

        return m_omega * std::sin(m_k_r * rr2 - m_omega * t);
    }
};

#endif /* INITIAL_CONDITONS_H_ */
