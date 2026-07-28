#ifndef RL_MATTER_PUMP_PARAMS_HPP_
#define RL_MATTER_PUMP_PARAMS_HPP_

#include "RLLumpState.hpp" // RL_MAX_LUMPS, RLLumpState

#include <array>
#include <cmath>

//! One pump "spotlight": a localized Gaussian driver centred on a lump.  The
//! centre is the lump's live tracked centroid (the spotlight follows the lump);
//! amplitude/frequency/phase are the RL agent's per-lump knobs.
struct RLPumpSite
{
    double center_x{0.0};
    double center_y{0.0};
    double center_z{0.0};
    double amplitude{0.0};
    double frequency{0.0};
    double phase{0.0};
    //! Sign of the field this spotlight drives.  +1 => canonical (Phi+) field,
    //! -1 => phantom (Phi-) field.  Used by the two-complex-field matter model
    //! to route the pump to the correct field; ignored by single-field models.
    int field_sign{1};
};

//! Multi-site pump: N independent spotlights sharing one envelope width and one
//! (global) L2_Ham safety governor.  Default-constructed => inactive.
struct RLMatterPumpParams
{
    RLPumpSite sites[RL_MAX_LUMPS]{};
    int num_sites{0};
    double width{1.5};
    double governor_center{0.035};
    double governor_width{0.003};
    double governor{1.0}; //!< host-precomputed tanh governor (device-safe)
    int num_fields{0};
    //! Closed-loop PD "trap" controller gains.  When ``k_p > 0`` the pump drives
    //! each field toward the target soliton at its moving trajectory centre with
    //! a proportional (k_p) + derivative (k_d) law instead of the legacy
    //! open-loop source.  The drive is error-proportional, hence self-limiting:
    //! it vanishes when the lump is on target, so it confines/transports the
    //! soliton without the over-driving that destroys structure.  ``k_p <= 0``
    //! falls back to the legacy additive-source pump.
    double k_p{0.0};
    double k_d{0.0};
    //! Per-sector gain override for the PHANTOM (Phi-) field.  Negative => the
    //! phantom sector inherits ``k_p`` / ``k_d`` (the historical behaviour, so
    //! existing configs are bit-identical).
    //!
    //! WHY THIS EXISTS.  The two sectors need different authority, not the same
    //! authority.  Measured on the mode-0 pump ladder (runs/pump_ladder_m0,
    //! t in [13, 30], confined_frac decay rate per unit time):
    //!
    //!     sector      pump off   pump on    outcome
    //!     canonical     0.100      ~0.000   decay ARRESTED (plateau ~0.50)
    //!     phantom       0.217       0.076   decay only SLOWED (2.9x)
    //!
    //! With one shared gain the controller removes a similar absolute amount of
    //! decay from each sector (~0.10 canonical, ~0.14 phantom), but the phantom
    //! sector's intrinsic dispersal rate is 2.2x larger, so the same authority
    //! arrests one and merely slows the other.  The shortfall is ~0.076 per unit
    //! time of unsuppressed phantom decay.  Raising only the phantom gain is the
    //! minimal lever that targets it without perturbing the canonical plateau.
    double k_p_phantom{-1.0};
    double k_d_phantom{-1.0};
    //! Trap-target matter profile shape and scale.  ``target_profile`` 0 =>
    //! Gaussian (legacy), 2 => sech bound lump (correct exponential tail).
    //! ``target_width`` is the physical bound-state size 1/sqrt(m^2-omega^2);
    //! when <= 0 the controller falls back to ``width``.  These must match the
    //! initial-data profile (PROFILE_SECH_BOUND, same width) so the controller
    //! drives the field toward the bound state it was initialised in, not a
    //! narrower (dispersing) Gaussian.
    int target_profile{0};
    double target_width{0.0};
    //! Trap-target central amplitude.  When > 0 the PD controller drives the
    //! field toward ``target_amp * envelope`` instead of ``site.amplitude *
    //! envelope``.  For a self-gravitating boson star this must be the star's
    //! central amplitude (bs_phi_c), NOT the transport knob ``well_depth`` --
    //! aiming at a much smaller amplitude makes the trap fight gravity and
    //! excites a breathing mode that collapses the star.  <= 0 => use
    //! ``site.amplitude`` (legacy).
    double target_amp{0.0};
};

namespace RLRuntime
{
inline double g_cached_l2_ham{0.0};

AMREX_FORCE_INLINE void publish_cached_L2_Ham(double l2_ham)
{
    g_cached_l2_ham = l2_ham;
}

AMREX_FORCE_INLINE double cached_L2_Ham() { return g_cached_l2_ham; }

//! Persistent, host-side lump state shared between the tracker (writer),
//! the observation collector (reader) and the pump-params builder (reader).
inline std::array<RLLumpState, RL_MAX_LUMPS> g_lump_state{};
inline int g_num_lumps{0};
inline bool g_lumps_seeded{false};

//! Seed the tracker's initial lump centres (once, before the evolution loop)
//! from the GRTresna elite geometry.
inline void seed_lumps(int num_lumps, const double *sx, const double *sy,
                       const double *sz)
{
    if (num_lumps > RL_MAX_LUMPS)
        num_lumps = RL_MAX_LUMPS;
    if (num_lumps < 0)
        num_lumps = 0;
    g_num_lumps = num_lumps;
    for (int k = 0; k < num_lumps; ++k)
    {
        g_lump_state[k]   = RLLumpState{};
        g_lump_state[k].x = sx[k];
        g_lump_state[k].y = sy[k];
        g_lump_state[k].z = sz[k];
    }
    g_lumps_seeded = true;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
tanh_governor(double cached_l2, double center, double width)
{
    const double margin = (center - cached_l2) / width;
    return 0.5 * (1.0 + std::tanh(margin));
}

//! Spatial × amplitude × governor weight of a single spotlight at (x,y,z).
//! No temporal factor (the caller applies cos/sin for the field component).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
compute_site_base(double x, double y, double z, const RLPumpSite &site,
                  double width, double governor)
{
    if (site.amplitude <= 0.0)
    {
        return 0.0;
    }
    const double dx = x - site.center_x;
    const double dy = y - site.center_y;
    const double dz = z - site.center_z;
    const double d2 = dx * dx + dy * dy + dz * dz;
    const double envelope = std::exp(-d2 / (2.0 * width * width));
    return site.amplitude * governor * envelope;
}

//! Pure Gaussian spotlight envelope at (x,y,z): exp(-r^2 / 2 width^2), with no
//! amplitude/governor factor.  Used by the closed-loop PD controller both as the
//! spatial localization window and to build the target soliton profile
//! (target |phi| = site.amplitude * envelope).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
compute_site_envelope(double x, double y, double z, const RLPumpSite &site,
                      double width, int profile = 0)
{
    const double dx = x - site.center_x;
    const double dy = y - site.center_y;
    const double dz = z - site.center_z;
    const double d2 = dx * dx + dy * dy + dz * dz;
    if (profile == 2)
    {
        // sech(r / width): bound-lump envelope matching the initial data
        // (PROFILE_SECH_BOUND).  Decays as 2 e^{-r/width} -- the correct tail,
        // so the controller does not erode the bound profile's wings the way a
        // Gaussian window (e^{-r^2}) would.
        const double r = std::sqrt(d2);
        return 1.0 / std::cosh(r / width);
    }
    return std::exp(-d2 / (2.0 * width * width));
}
} // namespace RLRuntime

#endif /* RL_MATTER_PUMP_PARAMS_HPP_ */
