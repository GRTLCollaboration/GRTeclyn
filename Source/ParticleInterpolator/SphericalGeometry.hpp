/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SPHERICALGEOMETRY_HPP_
#define SPHERICALGEOMETRY_HPP_

#include <AMReX.H>

// Other includes
#include <array>
#include <cmath>
#include <string>

//! This SurfaceGeometry template class provides spherical shell geometry
//! implementation for the SurfaceExtraction class
//! u = theta, v = phi
class SphericalGeometry
{
  private:
    std::array<amrex::ParticleReal, AMREX_SPACEDIM> m_center;

  public:
    SphericalGeometry(
        const std::array<amrex::ParticleReal, AMREX_SPACEDIM> &a_center)
        : m_center(a_center)
    {
    }

    //! returns the grid spacing in theta
    [[nodiscard]] static amrex::ParticleReal du(int a_num_points_theta)
    {
        return M_PI / (amrex::ParticleReal)(a_num_points_theta - 1);
    }

    //! returns the grid spacing in phi
    [[nodiscard]] static amrex::ParticleReal dv(int a_num_points_phi)
    {
        return 2.0 * M_PI / ((amrex::ParticleReal)a_num_points_phi);
    }

    //! returns the theta coordinate associated to the theta/u index
    [[nodiscard]] static amrex::ParticleReal u(int a_itheta,
                                               int a_num_points_theta)
    {
        return a_itheta * du(a_num_points_theta);
    }

    //! returns the phi coordinate associated to the phi/v index
    [[nodiscard]] static amrex::ParticleReal v(int a_iphi, int a_num_points_phi)
    {
        return a_iphi * dv(a_num_points_phi);
    }

    [[nodiscard]] static bool is_u_periodic() { return false; }
    [[nodiscard]] static bool is_v_periodic() { return true; }

    //! returns the Cartesian coordinate in direction a_dir with specified
    //! radius, theta and phi.
    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    [[nodiscard]] amrex::ParticleReal
    get_grid_coord(int a_dir, amrex::ParticleReal a_radius,
                   amrex::ParticleReal a_theta, amrex::ParticleReal a_phi) const
    // NOLINTEND(bugprone-easily-swappable-parameters)
    {
        switch (a_dir)
        {
        case (0):
            return m_center[0] + a_radius * sin(a_theta) * cos(a_phi);
        case (1):
            return m_center[1] + a_radius * sin(a_theta) * sin(a_phi);
        case (2):
            return m_center[2] + a_radius * cos(a_theta);
        default:
            amrex::Abort("SphericalGeometry: Direction not supported");
            return 0.;
        }
    }

    //! returns the area element on a sphere with radius a_radius at the point
    //! (a_theta, a_phi)
    [[nodiscard]] static amrex::ParticleReal
    area_element(amrex::ParticleReal a_radius, amrex::ParticleReal a_theta,
                 amrex::ParticleReal /*a_phi*/)
    {
        return a_radius * a_radius * sin(a_theta);
    }

    [[nodiscard]] static std::string param_name() { return "r"; }

    [[nodiscard]] static std::string u_name() { return "theta"; }

    [[nodiscard]] static std::string v_name() { return "phi"; }
};

#endif /* SPHERICALGEOMETRY_HPP_ */