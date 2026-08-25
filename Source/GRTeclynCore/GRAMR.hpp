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

// NOLINTNEXTLINE(cppcoreguidelines-special-member-functions)
class GRAMR : public amrex::Amr
{
    friend class GRAMRLevel;

  public:

    GRAMR(amrex::LevelBld *a_levelbld);
    ~GRAMR() override;

    void init(amrex::Real a_strt_time, amrex::Real a_stop_time) override;

    [[nodiscard]] double get_walltime_since_start() const;

    [[nodiscard]] double get_restart_time() const;

  private:

    void set_restart_time(double a_restart_time);

    double m_start_walltime{std::nan("0.0")};
    double m_restart_time{0.0};
};

#endif /* GRAMR_HPP_ */
