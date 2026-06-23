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
} // namespace RLRuntime

#endif /* RL_MATTER_PUMP_PARAMS_HPP_ */
