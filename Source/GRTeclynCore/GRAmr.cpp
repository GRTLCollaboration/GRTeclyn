/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#include "GRAmr.hpp"
#include "GRAmrLevel.hpp"
#include "SimulationParameters.hpp"

GRAmr::GRAmr(amrex::LevelBld *a_levelbld) : amrex::Amr(a_levelbld) {}

GRAmr::~GRAmr() = default;

void GRAmr::init(amrex::Real a_strt_time, amrex::Real a_stop_time)
{
    amrex::Amr::init(a_strt_time, a_stop_time);

    m_start_walltime = amrex::second();
}

amrex::Real GRAmr::get_walltime_since_start() const
{
    return amrex::second() - m_start_walltime;
}

amrex::Real GRAmr::get_restart_time() const { return m_restart_time; }

void GRAmr::set_restart_time(amrex::Real a_restart_time)
{
    m_restart_time = a_restart_time;
}