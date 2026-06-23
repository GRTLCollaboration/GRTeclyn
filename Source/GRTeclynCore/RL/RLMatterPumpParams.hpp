#ifndef RL_MATTER_PUMP_PARAMS_HPP_
#define RL_MATTER_PUMP_PARAMS_HPP_

#include <cmath>

//! Host/device pump parameters for Lump[0] RL forcing.
struct RLMatterPumpParams
{
    double amplitude{0.0};
    double frequency{0.0};
    double phase{0.0};
    double radius{5.0};
    double width{1.5};
    double governor_center{0.035};
    double governor_width{0.003};
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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
tanh_governor(double cached_l2, double center, double width)
{
    const double margin = (center - cached_l2) / width;
    return 0.5 * (1.0 + std::tanh(margin));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double compute_pump_drive(
    double x, double y, double z, double time,
    const RLMatterPumpParams &pump, double cached_l2)
{
    if (pump.num_fields < 1 || pump.amplitude <= 0.0)
    {
        return 0.0;
    }
    const double r = std::sqrt(x * x + y * y + z * z);
    const double envelope =
        std::exp(-std::pow(r - pump.radius, 2) /
                 (2.0 * pump.width * pump.width));
    const double governor =
        tanh_governor(cached_l2, pump.governor_center, pump.governor_width);
    return pump.amplitude * governor * envelope *
           std::cos(pump.frequency * time + pump.phase);
}
} // namespace RLRuntime

#endif /* RL_MATTER_PUMP_PARAMS_HPP_ */
