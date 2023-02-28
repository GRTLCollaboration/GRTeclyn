/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef GRAMR_HPP_
#define GRAMR_HPP_

// xxxxx#include "Lagrange.hpp"
#include "VariableType.hpp"
#include <AMReX_Amr.H>
#include <algorithm>
#include <chrono>
#include <ratio>
#include <vector>

/// A child of AMReX's AMR class to interface with tools which require
/// access to the whole AMR hierarchy
/**
 *It is necessary for many experimental features and allows us to
 *add said features later without breaking any user code.
 */

// Forward declaration for get_gramrlevels function declarations
class GRAMRLevel;

class SimulationParameters;

class GRAMR : public amrex::Amr
{
  private:
    using Clock = std::chrono::steady_clock;
    using Hours = std::chrono::duration<double, std::ratio<3600, 1>>;
    std::chrono::time_point<Clock> start_time = Clock::now();

  public:

    GRAMR(amrex::LevelBld *a_levelbld);
    virtual ~GRAMR();

    static void
    set_simulation_parameters(const SimulationParameters &a_sim_params);
    static SimulationParameters const &get_simulation_parameters();

  private:
    static SimulationParameters const *m_sim_params;

    // defined here due to auto return type
    auto get_walltime()
    {
        auto now      = Clock::now();
        auto duration = std::chrono::duration_cast<Hours>(now - start_time);

        return duration.count();
    }
};

#endif /* GRAMR_HPP_ */
