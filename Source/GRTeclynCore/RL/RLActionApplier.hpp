#ifndef RL_ACTION_APPLIER_HPP_
#define RL_ACTION_APPLIER_HPP_

#include "CCZ4RHS.hpp"

#include <algorithm>
#include <array>

//! Map normalized RL action targets to pump and gauge parameters (EMA).
inline void apply_rl_actions(
    double &pump_amplitude, double &pump_frequency, double &pump_phase,
    double pump_max_amplitude, CCZ4_params_t<> &ccz4_params,
    const std::array<double, 5> &actions)
{
    constexpr double k_ema = 0.2;

    const double requested_amp = std::clamp(actions[0], 0.0, pump_max_amplitude);
    pump_amplitude =
        (1.0 - k_ema) * pump_amplitude + k_ema * requested_amp;
    pump_amplitude = std::min(pump_amplitude, pump_max_amplitude);

    pump_frequency = std::clamp(actions[1], 0.0, 2.0);
    pump_phase     = actions[2];

    const double requested_lapse = std::clamp(1.0 + actions[3] * 0.5, 0.5, 1.5);
    const double requested_shift = std::clamp(0.75 + actions[4] * 0.25, 0.5, 1.0);

    ccz4_params.lapse_advec_coeff =
        (1.0 - k_ema) * ccz4_params.lapse_advec_coeff + k_ema * requested_lapse;
    ccz4_params.shift_Gamma_coeff =
        (1.0 - k_ema) * ccz4_params.shift_Gamma_coeff + k_ema * requested_shift;
}

#endif /* RL_ACTION_APPLIER_HPP_ */
