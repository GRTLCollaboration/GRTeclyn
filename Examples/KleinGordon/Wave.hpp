#ifndef WAVE_HPP_
#define WAVE_HPP_

// C++ std lib includes
#include <cmath>
// AMReX includes
#include <AMReX_MultiFab.H>
// GRTeclyn includes
#include "simd.hpp"
// KleinGordon includes
#include "StateVariables.hpp"

class Wave
{
  private:
    const amrex::Real m_k_r;
    const amrex::Real m_mass;
    const amrex::Real m_t0;

  public:
    Wave() = default;

    Wave(amrex::Real a_k_r, amrex::Real a_mass, amrex::Real a_initial_time)
        : m_k_r(a_k_r), m_mass(a_mass), m_t0(a_initial_time) {};

    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    calculate(const amrex::Real x, const amrex::Real y, const amrex::Real z,
              const amrex::Real t)
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
    calculate_derivative(const amrex::Real x, const amrex::Real y,
                         const amrex::Real z, const amrex::Real t)
    // NOLINTEND(bugprone-easily-swappable-parameters)
    {
        amrex::Real omega = m_k_r;

        // for the wave to be at the center of the grid, need to pass in
        // (x-x_midpt), (y-y_midpt) and (z-z_midpt)
        amrex::Real rr2 = x * x + y * y + z * z; // this is the radius

        return omega * std::sin(m_k_r * rr2 - omega * (t + m_t0));
    };

    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    template <class data_t, template <typename> class vars_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    compute_potential(data_t &V_of_phi, data_t &dVdphi,
                      const vars_t<data_t> &vars) const
    // NOLINTEND(bugprone-easily-swappable-parameters)
    {

        V_of_phi = 0.5 * m_mass * m_mass * vars.phi * vars.phi;

        dVdphi = m_mass * m_mass * vars.phi;
    }

    void calc_mf(amrex::MultiFab &mf_out, int dcomp,
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

                arrs[box_no](i, j, k, dcomp) = calculate(x, y, z, time);
                arrs[box_no](i, j, k, dcomp + 1) =
                    calculate_derivative(x, y, z, time);
            });
        amrex::Gpu::streamSynchronize();
    }
};
#endif // WAVE_HPP_
