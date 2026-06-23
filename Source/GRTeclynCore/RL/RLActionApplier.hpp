#ifndef RL_ACTION_APPLIER_HPP_
#define RL_ACTION_APPLIER_HPP_

#include "CCZ4RHS.hpp"
#include "GRTresnaScalarLayout.hpp" // GRTRESNA_MAX_INDEPENDENT_SCALARS

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

//! Map a normalized RL action vector to per-lump pump knobs + gauge (EMA).
//!
//! Action layout (length ``3 * num_lumps + 2``):
//!   [amp_0, freq_0, phase_0, amp_1, freq_1, phase_1, ..., lapse_advec, shift]
//! Pump amplitude uses EMA inertia + a hard cap; frequency/phase are set
//! directly (clamped).  Gauge coefficients also use EMA to avoid discontinuous
//! gauge jumps.
inline void apply_rl_actions(
    int num_lumps,
    std::array<double, GRTRESNA_MAX_INDEPENDENT_SCALARS> &pump_amplitude,
    std::array<double, GRTRESNA_MAX_INDEPENDENT_SCALARS> &pump_frequency,
    std::array<double, GRTRESNA_MAX_INDEPENDENT_SCALARS> &pump_phase,
    double pump_max_amplitude, CCZ4_params_t<> &ccz4_params,
    const std::vector<double> &actions)
{
    constexpr double k_ema = 0.2;

    if (num_lumps < 0)
        num_lumps = 0;
    if (num_lumps > GRTRESNA_MAX_INDEPENDENT_SCALARS)
        num_lumps = GRTRESNA_MAX_INDEPENDENT_SCALARS;

    const int needed = 3 * num_lumps + 2;
    if (static_cast<int>(actions.size()) < needed)
        return; // malformed action vector; leave state unchanged

    // All actions arrive raw-normalized in [-1, 1]; C++ is the single site that
    // maps them to physical bounds (and owns the hard safety caps).
    for (int s = 0; s < num_lumps; ++s)
    {
        const int base       = 3 * s;
        const double a_amp   = std::clamp(actions[base + 0], -1.0, 1.0);
        const double a_freq  = std::clamp(actions[base + 1], -1.0, 1.0);
        const double a_phase = std::clamp(actions[base + 2], -1.0, 1.0);
        const double requested_amp =
            0.5 * (a_amp + 1.0) * pump_max_amplitude; // [0, max]
        pump_amplitude[s] =
            (1.0 - k_ema) * pump_amplitude[s] + k_ema * requested_amp;
        pump_amplitude[s] = std::min(pump_amplitude[s], pump_max_amplitude);
        pump_frequency[s] = 0.5 * (a_freq + 1.0) * 2.0; // [0, 2]
        pump_phase[s]     = a_phase * M_PI;             // [-pi, pi]
    }

    const int g = 3 * num_lumps;
    const double requested_lapse =
        std::clamp(1.0 + actions[g + 0] * 0.5, 0.5, 1.5);
    const double requested_shift =
        std::clamp(0.75 + actions[g + 1] * 0.25, 0.5, 1.0);
    ccz4_params.lapse_advec_coeff =
        (1.0 - k_ema) * ccz4_params.lapse_advec_coeff + k_ema * requested_lapse;
    ccz4_params.shift_Gamma_coeff =
        (1.0 - k_ema) * ccz4_params.shift_Gamma_coeff + k_ema * requested_shift;
}

#endif /* RL_ACTION_APPLIER_HPP_ */
