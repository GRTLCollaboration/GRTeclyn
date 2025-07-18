#ifndef SINEGORDON_HPP_
#define SINEGORDON_HPP_

// C++ std lib includes
#include <cmath>
// AMReX includes
#include <AMReX_MultiFab.H>
// GRTeclyn includes
#include "simd.hpp"
// KleinGordon includes
#include "StateVariables.hpp"

class SineGordon
{
  private:
    const amrex::Real m_alpha;
    const amrex::Real m_t0;

  public:
    SineGordon() = default;

    SineGordon(amrex::Real a_alpha, amrex::Real a_initial_time)
        : m_alpha(a_alpha), m_t0(a_initial_time) {};

    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    calculate(const amrex::Real x, const amrex::Real t)
    // NOLINTEND(bugprone-easily-swappable-parameters)
    {
        // Sine Gordon 1D breather solution
        amrex::Real beta = std::sqrt(1.0 - m_alpha * m_alpha);

        return 4 * std::atan(beta * std::cos(m_alpha * (t + m_t0)) / m_alpha /
                             std::cosh(beta * x));
    };
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
    calculate_derivative(const amrex::Real x, const amrex::Real t)
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
              const amrex::Real t)
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

    void calc_mf_1d(amrex::MultiFab &mf_out, int dcomp,
                    const amrex::Geometry &geom, const amrex::Real time)
    {
        amrex::ParmParse pp;

        // Get the geometry of the problem
        const auto problo = geom.ProbLoArray();
        const auto dx     = geom.CellSizeArray();

        std::array<double, AMREX_SPACEDIM> center{};
        pp.query("center", center);

        auto const &arrs = mf_out.arrays();

        amrex::ParallelFor(
            mf_out, mf_out.nGrowVect(),
            [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
            {
                amrex::Real x = problo[0] + (i + 0.5) * dx[0] - center[0];

                arrs[box_no](i, j, k, dcomp) = calculate(x, time);
                arrs[box_no](i, j, k, dcomp + 1) =
                    calculate_derivative(x, time);
            });
        amrex::Gpu::streamSynchronize();
    }

    void calc_mf_3d(amrex::MultiFab &mf_out, int dcomp,
                    const amrex::Geometry &geom, const amrex::Real time)
    {
        amrex::ParmParse pp;

        // Get the geometry of the problem
        const auto problo = geom.ProbLoArray();
        const auto dx     = geom.CellSizeArray();

        std::array<double, AMREX_SPACEDIM> center{};
        pp.query("center", center);

        auto const &arrs = mf_out.arrays();

        amrex::ParallelFor(
            mf_out, mf_out.nGrowVect(),
            [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
            {
                amrex::Real x = problo[0] + (i + 0.5) * dx[0] - center[0];
                amrex::Real y = problo[1] + (j + 0.5) * dx[1] - center[1];
                amrex::Real z = problo[2] + (k + 0.5) * dx[2] - center[2];

                arrs[box_no](i, j, k, dcomp)     = calculate(x, y, z, time);
                arrs[box_no](i, j, k, dcomp + 1) = 0.0;
            });
        amrex::Gpu::streamSynchronize();
    }
};
#endif // SINEGORDON_HPP_
