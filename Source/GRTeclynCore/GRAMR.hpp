/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
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

// NOLINTNEXTLINE(cppcoreguidelines-special-member-functions)
class GRAMR : public amrex::Amr
{
    friend class GRAMRLevel;

  public:

    GRAMR(amrex::LevelBld *a_levelbld);
    ~GRAMR() override;

    virtual void init(amrex::Real a_strt_time,
                      amrex::Real a_stop_time) override;

    static void
    set_simulation_parameters(const SimulationParameters &a_sim_params);
    static const SimulationParameters &get_simulation_parameters();

    double get_walltime_since_start() const;

    double get_restart_time() const;

  private:
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
    static const SimulationParameters *m_sim_params;

    void set_restart_time(double a_restart_time);

    double m_start_walltime;
    double m_restart_time{0.0};
};

#endif /* GRAMR_HPP_ */
