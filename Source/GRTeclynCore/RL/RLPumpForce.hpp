#ifndef RL_PUMP_FORCE_HPP_
#define RL_PUMP_FORCE_HPP_

#include "RLMatterPumpParams.hpp"

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include <cmath>

//! Single source of truth for the PD / open-loop matter pump.
//!
//! SOLID: one responsibility — given pump params and the local field state,
//! compute the per-component sources added to ∂_t Π_A, and the equivalent
//! 3+1 force density (f_⊥, f_i) that those sources deposit into ∇_μ T^{μν}.
//!
//! Convention (matches the existing RHS and the pump_work diagnostic):
//!   rhs[c_Pi_A] += S_A
//!   f_⊥ = Σ_A s_A S_A Π_A
//!   f_i = Σ_A s_A S_A ∂_i φ_A
//! with gravitational sign s_A = +1 (canonical) or −1 (phantom).
struct RLPumpSources
{
    amrex::Real s1p{0.0}; //!< source on Π_1 of the + / single field
    amrex::Real s2p{0.0}; //!< source on Π_2 of the + / single field
    amrex::Real s1m{0.0}; //!< source on Π_1 of the − (phantom) field
    amrex::Real s2m{0.0}; //!< source on Π_2 of the − (phantom) field
};

struct RLPumpForceDensity
{
    amrex::Real f_perp{0.0};
    amrex::Real f_x{0.0};
    amrex::Real f_y{0.0};
    amrex::Real f_z{0.0};
};

namespace RLPumpForce
{
//! Accumulate PD / open-loop sources for sites matching `want_sign`.
//! `want_sign >= 0` selects canonical sites (field_sign >= 0);
//! `want_sign < 0` selects phantom sites. Pass `nullptr` for unused field
//! components (open-loop and PD then contribute 0 for those).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void accumulate_site_sources(
    amrex::Real &s1, amrex::Real &s2, const RLMatterPumpParams &pump,
    amrex::Real x, amrex::Real y, amrex::Real z, amrex::Real time,
    amrex::Real lapse, amrex::Real phi1, amrex::Real phi2, amrex::Real Pi1,
    amrex::Real Pi2, int want_sign)
{
    if (pump.num_sites < 1)
    {
        return;
    }
    const amrex::Real governor = pump.governor;
    // The PD-vs-open-loop MODE stays global (keyed on the shared ``k_p``), so a
    // zero phantom gain means "no PD drive on the phantom sector", never
    // "silently switch that sector to the legacy open-loop source".  Only the
    // gains inside the PD branch are per-sector.
    if (pump.k_p > 0.0)
    {
        const bool phantom   = (want_sign < 0);
        const amrex::Real kp = (phantom && pump.k_p_phantom >= 0.0)
                                   ? pump.k_p_phantom
                                   : pump.k_p;
        const amrex::Real kd = (phantom && pump.k_d_phantom >= 0.0)
                                   ? pump.k_d_phantom
                                   : pump.k_d;
        const amrex::Real inv_alpha = 1.0 / lapse;
        const amrex::Real tw =
            (pump.target_width > 0.0) ? pump.target_width : pump.width;
        // Superposed-target accumulators (see RLMatterPumpParams.hpp): the
        // sector's sites build ONE target field and one capped weight, and a
        // single PD error is taken at the end of the loop.
        amrex::Real env_sum{0.0};
        amrex::Real T1{0.0}, T2{0.0}, TP1{0.0}, TP2{0.0};
        for (int s = 0; s < pump.num_sites; ++s)
        {
            const auto &site = pump.sites[s];
            const bool is_plus = (site.field_sign >= 0);
            if (want_sign >= 0 ? !is_plus : is_plus)
            {
                continue;
            }
            if (site.amplitude <= 0.0)
            {
                continue;
            }
            const amrex::Real env = RLRuntime::compute_site_envelope(
                x, y, z, site, tw, pump.target_profile);
            if (env < 1.0e-8)
            {
                continue;
            }
            const amrex::Real amp_t =
                (pump.target_amp > 0.0) ? pump.target_amp : site.amplitude;
            const amrex::Real g     = amp_t * env;
            const amrex::Real arg   = -site.frequency * time + site.phase;
            const amrex::Real tphi1 = g * std::cos(arg);
            const amrex::Real tphi2 = g * std::sin(arg);
            const amrex::Real tPi1  = site.frequency * tphi2 * inv_alpha;
            const amrex::Real tPi2  = -site.frequency * tphi1 * inv_alpha;
            if (pump.superpose_targets != 0)
            {
                env_sum += env;
                T1 += tphi1;
                T2 += tphi2;
                TP1 += tPi1;
                TP2 += tPi2;
                continue;
            }
            const amrex::Real w = governor * env;
            s1 += w * (-kp * (phi1 - tphi1) - kd * (Pi1 - tPi1));
            s2 += w * (-kp * (phi2 - tphi2) - kd * (Pi2 - tPi2));
        }
        if (pump.superpose_targets != 0 && env_sum > 0.0)
        {
            const amrex::Real w =
                governor * ((env_sum < 1.0) ? env_sum : 1.0);
            s1 += w * (-kp * (phi1 - T1) - kd * (Pi1 - TP1));
            s2 += w * (-kp * (phi2 - T2) - kd * (Pi2 - TP2));
        }
        return;
    }

    // Legacy open-loop source pump (k_p <= 0).
    for (int s = 0; s < pump.num_sites; ++s)
    {
        const auto &site = pump.sites[s];
        const bool is_plus = (site.field_sign >= 0);
        if (want_sign >= 0 ? !is_plus : is_plus)
        {
            continue;
        }
        const amrex::Real base =
            RLRuntime::compute_site_base(x, y, z, site, pump.width, governor);
        if (base <= 0.0)
        {
            continue;
        }
        const amrex::Real arg = site.frequency * time + site.phase;
        s1 += base * std::cos(arg);
        s2 += base * std::sin(arg);
    }
}

//! Sources for a single complex field (all sites drive the + slot; field_sign
//! is ignored for routing — matches ComplexScalarField / ComplexExotic).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE RLPumpSources
compute_single_field_sources(const RLMatterPumpParams &pump, amrex::Real x,
                             amrex::Real y, amrex::Real z, amrex::Real time,
                             amrex::Real lapse, amrex::Real phi1,
                             amrex::Real phi2, amrex::Real Pi1, amrex::Real Pi2)
{
    RLPumpSources out;
    if (pump.num_sites < 1)
    {
        return out;
    }
    // Temporarily treat every site as canonical for routing.
    // Sites keep their stored field_sign; we accept both by calling twice
    // would double-count, so we inline a sign-agnostic loop here.
    const amrex::Real governor = pump.governor;
    if (pump.k_p > 0.0)
    {
        const amrex::Real kp        = pump.k_p;
        const amrex::Real kd        = pump.k_d;
        const amrex::Real inv_alpha = 1.0 / lapse;
        const amrex::Real tw =
            (pump.target_width > 0.0) ? pump.target_width : pump.width;
        for (int s = 0; s < pump.num_sites; ++s)
        {
            const auto &site = pump.sites[s];
            if (site.amplitude <= 0.0)
            {
                continue;
            }
            const amrex::Real env = RLRuntime::compute_site_envelope(
                x, y, z, site, tw, pump.target_profile);
            if (env < 1.0e-8)
            {
                continue;
            }
            const amrex::Real amp_t =
                (pump.target_amp > 0.0) ? pump.target_amp : site.amplitude;
            const amrex::Real g     = amp_t * env;
            const amrex::Real arg   = -site.frequency * time + site.phase;
            const amrex::Real tphi1 = g * std::cos(arg);
            const amrex::Real tphi2 = g * std::sin(arg);
            const amrex::Real tPi1  = site.frequency * tphi2 * inv_alpha;
            const amrex::Real tPi2  = -site.frequency * tphi1 * inv_alpha;
            const amrex::Real w     = governor * env;
            out.s1p += w * (-kp * (phi1 - tphi1) - kd * (Pi1 - tPi1));
            out.s2p += w * (-kp * (phi2 - tphi2) - kd * (Pi2 - tPi2));
        }
        return out;
    }
    for (int s = 0; s < pump.num_sites; ++s)
    {
        const amrex::Real base = RLRuntime::compute_site_base(
            x, y, z, pump.sites[s], pump.width, governor);
        if (base <= 0.0)
        {
            continue;
        }
        const amrex::Real arg =
            pump.sites[s].frequency * time + pump.sites[s].phase;
        out.s1p += base * std::cos(arg);
        out.s2p += base * std::sin(arg);
    }
    return out;
}

//! Sources for the two-complex-field model (routes by site.field_sign).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE RLPumpSources
compute_bicomplex_sources(const RLMatterPumpParams &pump, amrex::Real x,
                          amrex::Real y, amrex::Real z, amrex::Real time,
                          amrex::Real lapse, amrex::Real phi1p,
                          amrex::Real phi2p, amrex::Real Pi1p, amrex::Real Pi2p,
                          amrex::Real phi1m, amrex::Real phi2m, amrex::Real Pi1m,
                          amrex::Real Pi2m)
{
    // accumulate_site_sources filters by site.field_sign internally in BOTH
    // the PD and the open-loop branch, so one pass per sign is exact and
    // double-counts nothing.
    RLPumpSources out;
    accumulate_site_sources(out.s1p, out.s2p, pump, x, y, z, time, lapse, phi1p,
                            phi2p, Pi1p, Pi2p, /*want_sign=*/+1);
    accumulate_site_sources(out.s1m, out.s2m, pump, x, y, z, time, lapse, phi1m,
                            phi2m, Pi1m, Pi2m, /*want_sign=*/-1);
    return out;
}

//! 3+1 force density from sources and local field / derivatives.
//! `sign_p` / `sign_m` are the gravitational couplings (+1 / −1).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE RLPumpForceDensity
force_density_from_sources(const RLPumpSources &src, amrex::Real sign_p,
                           amrex::Real Pi1p, amrex::Real Pi2p, amrex::Real dphi1p_x,
                           amrex::Real dphi1p_y, amrex::Real dphi1p_z,
                           amrex::Real dphi2p_x, amrex::Real dphi2p_y,
                           amrex::Real dphi2p_z, amrex::Real sign_m,
                           amrex::Real Pi1m, amrex::Real Pi2m, amrex::Real dphi1m_x,
                           amrex::Real dphi1m_y, amrex::Real dphi1m_z,
                           amrex::Real dphi2m_x, amrex::Real dphi2m_y,
                           amrex::Real dphi2m_z)
{
    RLPumpForceDensity f;
    f.f_perp += sign_p * (src.s1p * Pi1p + src.s2p * Pi2p);
    f.f_x += sign_p * (src.s1p * dphi1p_x + src.s2p * dphi2p_x);
    f.f_y += sign_p * (src.s1p * dphi1p_y + src.s2p * dphi2p_y);
    f.f_z += sign_p * (src.s1p * dphi1p_z + src.s2p * dphi2p_z);
    f.f_perp += sign_m * (src.s1m * Pi1m + src.s2m * Pi2m);
    f.f_x += sign_m * (src.s1m * dphi1m_x + src.s2m * dphi2m_x);
    f.f_y += sign_m * (src.s1m * dphi1m_y + src.s2m * dphi2m_y);
    f.f_z += sign_m * (src.s1m * dphi1m_z + src.s2m * dphi2m_z);
    return f;
}

//! Convenience: single-field force density (phantom slot unused).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE RLPumpForceDensity
force_density_single(const RLPumpSources &src, amrex::Real sign, amrex::Real Pi1,
                     amrex::Real Pi2, amrex::Real dphi1_x, amrex::Real dphi1_y,
                     amrex::Real dphi1_z, amrex::Real dphi2_x,
                     amrex::Real dphi2_y, amrex::Real dphi2_z)
{
    return force_density_from_sources(src, sign, Pi1, Pi2, dphi1_x, dphi1_y,
                                      dphi1_z, dphi2_x, dphi2_y, dphi2_z,
                                      /*sign_m=*/0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                      0.0, 0.0, 0.0);
}
} // namespace RLPumpForce

#endif /* RL_PUMP_FORCE_HPP_ */
