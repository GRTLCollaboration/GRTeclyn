/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef BHMOVINGPUNCTUREGAUGE_HPP_
#define BHMOVINGPUNCTUREGAUGE_HPP_

#include "MovingPunctureGauge.hpp"

#include <AMReX_BLassert.H>

#include <array>

/// A moving puncture gauge with a mass-dependent eta for binary black holes.
template <class deriv_t = FourthOrderDerivatives>
class BHMovingPunctureGauge : public MovingPunctureGauge<deriv_t>
{
  public:
    static constexpr std::size_t num_punctures = 2;
    static constexpr std::size_t num_puncture_coords =
        num_punctures * AMREX_SPACEDIM;
    using puncture_masses_t = std::array<amrex::Real, num_punctures>;
    using puncture_coords_t = std::array<amrex::Real, num_puncture_coords>;

  private:
    puncture_masses_t m_puncture_masses{};
    puncture_coords_t m_puncture_coords{};

  public:
    template <class masses_t, class coords_t>
    BHMovingPunctureGauge(amrex::Real a_dx, const masses_t &a_puncture_masses,
                          const coords_t &a_puncture_coords)
        : MovingPunctureGauge<deriv_t>(a_dx)
    {
        for (std::size_t puncture = 0; puncture < num_punctures; ++puncture)
        {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                a_puncture_masses[puncture] > 0.0,
                "Puncture masses must be greater than zero when calculating "
                "the puncture-dependent gauge eta");
            m_puncture_masses[puncture] = a_puncture_masses[puncture];
        }
        for (std::size_t coord = 0; coord < num_puncture_coords; ++coord)
        {
            m_puncture_coords[coord] = a_puncture_coords[coord];
        }
    }

    /// Compute the spatially varying Gamma-driver damping parameter.
    /** First calculates \f$\eta_p=1/(2m_p)\f$ and interpolates between the
     * two punctures using
     * \f$\eta_*=(\eta_1r_2^2+\eta_2r_1^2)/(r_1^2+r_2^2)\f$. It then applies
     * the optional far-field cutoff configured by the inherited gauge
     * parameters.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
    compute_eta(amrex::Real &eta_of_x, int ix, int iy, int iz) const
    {
        const amrex::IntVect cell_index(ix, iy, iz);
        const Coordinates coords(cell_index, this->m_dx, this->m_params.center);
        const amrex::Real radius = coords.get_radius();
        const amrex::Real eta_cutoff_radius_squared =
            this->m_params.eta_cutoff_radius * this->m_params.eta_cutoff_radius;
        const amrex::Real eta_cutoff_factor =
            eta_cutoff_radius_squared /
            (radius * radius + eta_cutoff_radius_squared);

        amrex::Real radius_squared[num_punctures]{};
        for (std::size_t puncture = 0; puncture < num_punctures; ++puncture)
        {
            const std::size_t first_coord = puncture * AMREX_SPACEDIM;
            const std::array<amrex::Real, AMREX_SPACEDIM> puncture_center{
                m_puncture_coords[first_coord],
                m_puncture_coords[first_coord + 1],
                m_puncture_coords[first_coord + 2]};
            const Coordinates puncture_coords(cell_index, this->m_dx,
                                              puncture_center);
            radius_squared[puncture] = puncture_coords.x * puncture_coords.x +
                                       puncture_coords.y * puncture_coords.y +
                                       puncture_coords.z * puncture_coords.z;
        }

        const amrex::Real radius_squared_sum =
            radius_squared[0] + radius_squared[1];
        amrex::Real central_eta =
            0.25 * (1.0 / m_puncture_masses[0] + 1.0 / m_puncture_masses[1]);
        if (radius_squared_sum > 0.0)
        {
            central_eta = 0.5 *
                          (radius_squared[1] / m_puncture_masses[0] +
                           radius_squared[0] / m_puncture_masses[1]) /
                          radius_squared_sum;
        }

        eta_of_x =
            central_eta * (this->m_params.eta_cutoff_coeff * eta_cutoff_factor +
                           (1.0 - this->m_params.eta_cutoff_coeff));
    }

    /// Calculate the gauge RHS using the mass-dependent eta profile.
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    calculate_rhs(int ix, int iy, int iz, const amrex::Array4<amrex::Real> &rhs,
                  const amrex::Array4<const amrex::Real> &state) const
    {
        amrex::Real eta_of_x{};
        compute_eta(eta_of_x, ix, iy, iz);
        this->calculate_rhs_with_eta(ix, iy, iz, rhs, state, eta_of_x);
    }
};

#endif /* BHMOVINGPUNCTUREGAUGE_HPP_ */
