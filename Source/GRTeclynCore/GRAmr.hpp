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
#include <limits>
#include <ratio>
#include <vector>

/// A child of AMReX's AMR class to interface with tools which require
/// access to the whole AMR hierarchy
/**
 *It is necessary for many experimental features and allows us to
 *add said features later without breaking any user code.
 */

// Forward declaration for get_gramrlevels function declarations
class GRAmrLevel;

// NOLINTNEXTLINE(cppcoreguidelines-special-member-functions)
class GRAmr : public amrex::Amr
{
    friend class GRAmrLevel;

  public:

    GRAmr(amrex::LevelBld *a_levelbld);
    ~GRAmr() override;

    void init(amrex::Real a_strt_time, amrex::Real a_stop_time) override;

    [[nodiscard]] amrex::Real get_walltime_since_start() const;

    [[nodiscard]] amrex::Real get_restart_time() const;

  private:

    void set_restart_time(amrex::Real a_restart_time);

    // NaN marks the start time as unset and matches the configured Real type.
    amrex::Real m_start_walltime{std::numeric_limits<amrex::Real>::quiet_NaN()};
    amrex::Real m_restart_time{0.0};
};

#endif /* GRAMR_HPP_ */
