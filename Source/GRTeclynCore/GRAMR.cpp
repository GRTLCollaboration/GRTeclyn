/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#include "GRAMR.hpp"
#include "GRAMRLevel.hpp"
#include "SimulationParameters.hpp"

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
const SimulationParameters *GRAMR::m_sim_params = nullptr;

GRAMR::GRAMR(amrex::LevelBld *a_levelbld) : amrex::Amr(a_levelbld) {}

GRAMR::~GRAMR() = default;

void GRAMR::set_simulation_parameters(const SimulationParameters &a_sim_params)
{
    m_sim_params = &a_sim_params;
}

const SimulationParameters &GRAMR::get_simulation_parameters()
{
    return *m_sim_params;
}

void GRAMR::init(amrex::Real a_strt_time, amrex::Real a_stop_time)
{
    amrex::Amr::init(a_strt_time, a_stop_time);

    m_start_walltime = amrex::second();
}

double GRAMR::get_walltime_since_start() const
{
    return amrex::second() - m_start_walltime;
}

double GRAMR::get_restart_time() const { return m_restart_time; }

void GRAMR::set_restart_time(double a_restart_time)
{
    m_restart_time = a_restart_time;
}

// we do not need std::unique_ptr here strictly speaking. We can actually omit
// this step.
void GRAMR::convert_derived_multifabs(
    const amrex::Vector<std::unique_ptr<amrex::MultiFab>> &inputs,
    amrex::Vector<amrex::MultiFab *> &fields)
{
    fields.clear();
    fields.reserve(inputs.size());
    for (auto const &level_content : inputs)
    {
        fields.push_back(level_content.get());
    }
}