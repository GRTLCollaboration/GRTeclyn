/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef WAVE_HPP_
#define WAVE_HPP_

// C++ std lib includes
#include <cmath>
// AMReX includes
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
// KleinGordon includes
#include "GRParmParse.hpp"
#include "StateVariables.hpp"

class Wave
{
  public:

    struct params_t
    {
        amrex::Real k_r{1.0};
        amrex::Real scalar_mass{0.0};

        static void check_params()
        {
            GRParmParse wave_pp("wave");
            amrex::Real k_r{1.0};
            wave_pp.queryAdd("wave_vector", k_r);
            amrex::Real scalar_mass{0.0};
            wave_pp.queryAdd(
                "scalar_mass",
                scalar_mass); // What is the mass of the scalar particle?
        }

        void fill_params()
        {
            GRParmParse wave_pp("wave");
            wave_pp.get("wave_vector", k_r);
            wave_pp.get("scalar_mass", scalar_mass);
        }
    };

    params_t m_params;
    amrex::Real m_t0{0.0};

    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
    Wave()
    {
        amrex::ParmParse pp;
        pp.query("klein_gordon.initial_time", m_t0);

        m_params.fill_params();
    }

    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    calculate(const amrex::Real x, const amrex::Real y, const amrex::Real z,
              const amrex::Real t) const
    // NOLINTEND(bugprone-easily-swappable-parameters)
    {
        amrex::Real omega = m_params.k_r;

        // for the wave to be at the center of the grid, need to pass in
        // (x-x_midpt), (y-y_midpt) and (z-z_midpt)
        amrex::Real rr2 = x * x + y * y + z * z; // this is the radius

        return std::cos(m_params.k_r * rr2 - omega * (t + m_t0));
    };

    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    calculate_time_derivative(const amrex::Real x, const amrex::Real y,
                              const amrex::Real z, const amrex::Real t) const
    // NOLINTEND(bugprone-easily-swappable-parameters)
    {
        amrex::Real omega = m_params.k_r;

        // for the wave to be at the center of the grid, need to pass in
        // (x-x_midpt), (y-y_midpt) and (z-z_midpt)
        amrex::Real rr2 = x * x + y * y + z * z; // this is the radius

        return omega * std::sin(m_params.k_r * rr2 - omega * (t + m_t0));
    };

    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    compute_potential(amrex::Real &V_of_phi, amrex::Real &dVdphi,
                      const amrex::Real &phi) const
    // NOLINTEND(bugprone-easily-swappable-parameters)
    {

        V_of_phi =
            0.5 * m_params.scalar_mass * m_params.scalar_mass * phi * phi;

        dVdphi = m_params.scalar_mass * m_params.scalar_mass * phi;
    }
};

#endif // WAVE_HPP_
