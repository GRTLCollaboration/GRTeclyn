/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef POTENTIAL_HPP_
#define POTENTIAL_HPP_

#include "GRParmParse.hpp"
#include "ScalarFieldVars.hpp"

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

class Potential
{
  public:
    struct params_t
    {
        amrex::Real scalar_mass{1.0};

        static void check_params()
        {
            GRParmParse scalar_field_pp("scalar_field");
            amrex::Real scalar_mass{1.0};
            scalar_field_pp.queryAdd("scalar_mass", scalar_mass);
            if (scalar_mass < 0.0)
            {
                scalar_field_pp.error("scalar_mass", "must be >= 0.0");
            }

            GRParmParse geometry_pp("geometry");
            amrex::Real coarsest_dx{};
            geometry_pp.get("coarsest_dx", coarsest_dx);

            GRParmParse evolution_pp("evolution");
            amrex::Real dt_multiplier{};
            evolution_pp.get("dt_multiplier", dt_multiplier);
            if (scalar_mass >= 0.2 / coarsest_dx / dt_multiplier)
            {
                scalar_field_pp.warning(
                    "scalar_mass",
                    "oscillations of the scalar field may not be resolved on "
                    "the coarsest level");
            }
        }

        void fill_params()
        {
            GRParmParse scalar_field_pp("scalar_field");
            scalar_field_pp.get("scalar_mass", scalar_mass);
        }
    };

    Potential() { m_params.fill_params(); }

    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE explicit Potential(params_t a_params)
        : m_params(a_params)
    {
    }

    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    compute_potential(amrex::Real &V_of_phi, amrex::Real &dVdphi,
                      const ScalarFieldVars &vars) const
    {
        const amrex::Real mass_times_phi = m_params.scalar_mass * vars.phi();
        V_of_phi = 0.5 * mass_times_phi * mass_times_phi;
        dVdphi   = m_params.scalar_mass * m_params.scalar_mass * vars.phi();
    }
    // NOLINTEND(bugprone-easily-swappable-parameters)

  private:
    params_t m_params{};
};

#endif /* POTENTIAL_HPP_ */
