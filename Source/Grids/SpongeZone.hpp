#ifndef SPONGEZONE_HPP_
#define SPONGEZONE_HPP_

#include "FourthOrderDerivatives.hpp"
#include <AMReX_Array4.H>
#include <AMReX_REAL.H>
#include <array>
#include <cmath>

/// Radially-ramped numerical sponge zone.
///
/// Adds spatially varying extra Kreiss-Oliger dissipation in a spherical shell
/// between `inner_radius` and `outer_radius` around `center`.  The ramp
/// function is a quartic (power=4) profile by default:
///
///   f(r) = 0                                        for r <= inner_radius
///   f(r) = ((r - inner) / (outer - inner))^power    for inner < r < outer
///   f(r) = 1                                        for r >= outer_radius
///
///   extra_sigma(r) = sponge_strength * f(r)
///
/// The sponge is applied as an ADDITIONAL dissipation pass after the main
/// CCZ4RHSWithMatter evaluation, so the total effective dissipation is:
///   sigma_total(r) = sigma_base + extra_sigma(r)
///
/// Design: Option A from SPONGE_ZONE_INTEGRATION_PLAN.md — per-example,
/// applied in RadialRecipeMatterDispatch::eval_rhs() as a second ParallelFor.
/// No changes to shared CCZ4 infrastructure.
struct SpongeZoneParams
{
    bool enabled{false};
    double inner_radius{24.0};
    double outer_radius{32.0};
    double strength{4.0}; //!< Extra sigma at full ramp (added on top of base sigma)
    int ramp_power{4};    //!< Ramp exponent (4 = quartic)
    std::array<double, AMREX_SPACEDIM> center{};
};

/// GPU-capturable sponge operator.  Constructed once per RHS call from the
/// params and dx, then applied via ParallelFor.
struct SpongeZone
{
    bool enabled;
    double inner_radius;
    double inv_width; //!< 1.0 / (outer - inner)
    double strength;
    int ramp_power;
    double dx;
    std::array<double, AMREX_SPACEDIM> center;

    SpongeZone() = default;

    explicit SpongeZone(const SpongeZoneParams &p, double a_dx)
        : enabled(p.enabled), inner_radius(p.inner_radius),
          inv_width((p.outer_radius > p.inner_radius)
                        ? 1.0 / (p.outer_radius - p.inner_radius)
                        : 0.0),
          strength(p.strength), ramp_power(p.ramp_power), dx(a_dx),
          center(p.center)
    {
    }

    /// Apply extra dissipation at point (i,j,k).  Adds to existing RHS.
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    apply(int i, int j, int k,
          const amrex::Array4<amrex::Real> &rhs,
          const amrex::Array4<amrex::Real const> &state) const
    {
        if (!enabled)
            return;

        // Compute physical position relative to center.
        double x = (static_cast<double>(i) + 0.5) * dx - center[0];
        double y = (static_cast<double>(j) + 0.5) * dx - center[1];
        double z = (static_cast<double>(k) + 0.5) * dx - center[2];
        double r = std::sqrt(x * x + y * y + z * z);

        if (r <= inner_radius)
            return;

        // Ramp function: ((r - inner) / width)^power, capped at 1.
        double s = (r - inner_radius) * inv_width;
        if (s > 1.0)
            s = 1.0;
        // General power ramp via std::pow.
        double f = std::pow(s, ramp_power);

        double extra_sigma = strength * f;

        // Apply dissipation stencil for each variable via the same API used
        // by CCZ4RHSWithMatter: Array4::cellData() → CellData reference.
        FourthOrderDerivatives deriv(dx);
        const amrex::CellData<amrex::Real> &rhs_cell = rhs.cellData(i, j, k);
        deriv.add_dissipation(i, j, k, rhs_cell, state, extra_sigma);
    }
};

#endif /* SPONGEZONE_HPP_ */
