/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
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
