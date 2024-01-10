/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#include "GRAMRLevel.hpp"

void GRAMRLevel::variableCleanUp()
{
    desc_lst.clear();
    derive_lst.clear();
}

GRAMRLevel::GRAMRLevel() = default;

GRAMRLevel::GRAMRLevel(amrex::Amr &papa, int lev, const amrex::Geometry &geom,
                       const amrex::BoxArray &box_array,
                       const amrex::DistributionMapping &distribution_mapping,
                       amrex::Real time)
    : amrex::AmrLevel(papa, lev, geom, box_array, distribution_mapping, time),
      m_num_ghosts(simParams().num_ghosts)
{

    m_boundaries.define(simParams().center, simParams().boundary_params, geom,
                        m_num_ghosts);
}

GRAMRLevel::~GRAMRLevel() = default;

const SimulationParameters &GRAMRLevel::simParams()
{
    return GRAMR::get_simulation_parameters();
}

GRAMR *GRAMRLevel::get_gramr_ptr()
{
    if (m_gramr_ptr == nullptr)
    {
        if (parent == nullptr)
        {
            amrex::Abort("AmrLevel::parent is null");
        }
        m_gramr_ptr = dynamic_cast<GRAMR *>(parent);
    }
    return m_gramr_ptr;
}

void GRAMRLevel::computeInitialDt(
    int finest_level, int /*sub_cycle*/, amrex::Vector<int> & /*n_cycle*/,
    const amrex::Vector<amrex::IntVect> & /*ref_ratio*/,
    amrex::Vector<amrex::Real> &dt_level, amrex::Real /*stop_time*/)
{
    // Level 0 will do it for all levels
    if (Level() == 0)
    {
        double dt_multiplier = simParams().dt_multiplier;
        for (int i = 0; i <= finest_level; ++i)
        {
            dt_level[i] = dt_multiplier * parent->Geom(i).CellSize(0);
        }
    }
}

void GRAMRLevel::computeNewDt(
    int finest_level, int /*sub_cycle*/, amrex::Vector<int> & /*n_cycle*/,
    const amrex::Vector<amrex::IntVect> & /*ref_ratio*/,
    amrex::Vector<amrex::Real> &dt_min, amrex::Vector<amrex::Real> &dt_level,
    amrex::Real /*stop_time*/, int /*post_regrid_flag*/)
{
    // This is called at the end of a coarse time step
    // Level 0 will do it for all levels
    if (Level() == 0)
    {
        double dt_multiplier = simParams().dt_multiplier;
        for (int i = 0; i <= finest_level; ++i)
        {
            dt_min[i] = dt_level[i] =
                dt_multiplier * parent->Geom(i).CellSize(0);
        }
    }
}

amrex::Real GRAMRLevel::advance(amrex::Real time, amrex::Real dt, int iteration,
                                int ncycle)
{
    BL_PROFILE("GRAMRLevel::advance()");
    double seconds_per_hour = 3600;
    double evolution_speed  = (time - get_gramr_ptr()->get_restart_time()) *
                             seconds_per_hour /
                             get_gramr_ptr()->get_walltime_since_start();
    amrex::Print() << "[Level " << Level() << " step "
                   << parent->levelSteps(Level()) + 1
                   << "] average evolution speed = " << evolution_speed
                   << " code units/h\n";

    for (int k = 0; k < NUM_STATE_TYPE; k++)
    {
        state[k].allocOldData();
        state[k].swapTimeLevels(dt);
    }

    amrex::AmrLevel::RK(
        4, State_Type, time, dt, iteration, ncycle,
        [&](int /*stage*/, amrex::MultiFab &rhs, const amrex::MultiFab &soln,
            amrex::Real t, amrex::Real /*dtsub*/)
        {
            // NOLINTNEXTLINE(cppcoreguidelines-pro-type-const-cast)
            specificEvalRHS(const_cast<amrex::MultiFab &>(soln), rhs, t);
            m_boundaries.apply_sommerfeld_boundaries(rhs, soln);
        },
        [&](int /*stage*/, amrex::MultiFab &soln) { specificUpdateODE(soln); });

    specificAdvance();

    return dt;
}

void GRAMRLevel::post_timestep(int /*iteration*/)
{
    BL_PROFILE("GRAMRLevel::post_timestep()");
    const int lev = Level();
    if (lev < parent->finestLevel())
    {
        auto &fine_level        = parent->getLevel(Level() + 1);
        amrex::MultiFab &S_fine = fine_level.get_new_data(State_Type);
        amrex::MultiFab &S_crse = this->get_new_data(State_Type);
        amrex::Real t           = get_state_data(State_Type).curTime();

        amrex::IntVect ratio = parent->refRatio(lev);
        AMREX_ASSERT(ratio == 2 || ratio == 4);
        if (ratio == 2)
        {
            // Need to fill one ghost cell for the high-order interpolation
            // below
            FillPatch(fine_level, S_fine, 1, t, State_Type, 0, S_fine.nComp());
        }

        FourthOrderInterpFromFineToCoarse(S_crse, 0, 2, S_fine, ratio);
    }

    specificPostTimeStep();
}

void GRAMRLevel::post_regrid(int /*lbase*/, int /*new_finest*/)
{
    // xxxxx Do we need to do anything after regrid?
}

void GRAMRLevel::post_init(amrex::Real /*stop_time*/)
{
    if (Level() == 0)
    {
        get_gramr_ptr()->set_restart_time(get_gramr_ptr()->cumTime());
    }
}

void GRAMRLevel::post_restart()
{
    if (Level() == 0)
    {
        get_gramr_ptr()->set_restart_time(get_gramr_ptr()->cumTime());
    }
}

void GRAMRLevel::init(amrex::AmrLevel &old)
{
    BL_PROFILE("GRAMRLevel::init()");
    amrex::Real dt_new    = parent->dtLevel(level);
    amrex::Real cur_time  = old.get_state_data(State_Type).curTime();
    amrex::Real prev_time = old.get_state_data(State_Type).prevTime();
    amrex::Real dt_old    = cur_time - prev_time;
    setTimeLevel(cur_time, dt_old, dt_new);

    amrex::MultiFab &S_new = get_new_data(State_Type);
    FillPatch(old, S_new, 0, cur_time, State_Type, 0, S_new.nComp());
}

void GRAMRLevel::init()
{
    BL_PROFILE("GRAMRLevel::init()");
    amrex::Real dt = parent->dtLevel(level);
    const auto &coarse_state =
        parent->getLevel(level - 1).get_state_data(State_Type);
    amrex::Real cur_time  = coarse_state.curTime();
    amrex::Real prev_time = coarse_state.prevTime();
    amrex::Real dt_old =
        (cur_time - prev_time) /
        static_cast<amrex::Real>(parent->MaxRefRatio(level - 1));
    setTimeLevel(cur_time, dt_old, dt);

    amrex::MultiFab &S_new = get_new_data(State_Type);
    FillCoarsePatch(S_new, 0, cur_time, State_Type, 0, S_new.nComp());
}

void GRAMRLevel::writePlotFilePre(const std::string & /*dir*/,
                                  std::ostream & /*os*/)
{
    m_is_writing_plotfile = true;
}

void GRAMRLevel::writePlotFilePost(const std::string & /*dir*/,
                                   std::ostream & /*os*/)
{
    m_is_writing_plotfile = false;
}

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
std::unique_ptr<amrex::MultiFab> GRAMRLevel::derive(const std::string &name,
                                                    amrex::Real time, int ngrow)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    std::unique_ptr<amrex::MultiFab> multifab;
    const amrex::DeriveRec *rec = derive_lst.get(name);
    if (rec != nullptr)
    {
        multifab = std::make_unique<amrex::MultiFab>(
            this->boxArray(), this->DistributionMap(), rec->numState(), ngrow);
        derive(name, time, *multifab, 0);
    }
    else
    {
        amrex::Abort("Unknown derived variable");
    }
    return multifab;
}

void GRAMRLevel::derive(const std::string &name, amrex::Real time,
                        amrex::MultiFab &multifab, int dcomp)
{
    amrex::Abort("GRAMRLevel::derive(): Implement this function in the child "
                 "level class if derived variables are required.");
}