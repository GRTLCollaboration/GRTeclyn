/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef POTENTIAL_HPP_
#define POTENTIAL_HPP_

#include "ScalarFieldVars.hpp"

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

class Potential
{
  public:
    struct params_t
    {
        amrex::Real scalar_mass{1.0};
    };

    Potential() = default;

    AMREX_GPU_HOST_DEVICE
        AMREX_FORCE_INLINE explicit Potential(params_t a_params)
        : m_params(a_params)
    {
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    compute_potential(amrex::Real &V_of_phi, amrex::Real &dVdphi,
                      const ScalarFieldVars &vars) const
    {
        const amrex::Real mass_times_phi = m_params.scalar_mass * vars.phi();
        V_of_phi = 0.5 * mass_times_phi * mass_times_phi;
        dVdphi   = m_params.scalar_mass * m_params.scalar_mass * vars.phi();
    }

  private:
    params_t m_params{};
};

#endif /* POTENTIAL_HPP_ */
