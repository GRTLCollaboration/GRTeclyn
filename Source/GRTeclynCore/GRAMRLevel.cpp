/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#include "GRAMRLevel.hpp"
#include "NullBCFill.hpp"

void GRAMRLevel::stateVariableSetUp()
{
    const int nghost = simParams().num_ghosts;
    desc_lst.addDescriptor(zero_state_index, amrex::IndexType::TheCellType(),
                           amrex::StateDescriptor::Point, nghost, NUM_VARS,
                           &amrex::cell_quartic_interp);

    BoundaryConditions::params_t bparms = simParams().boundary_params;
    BoundaryConditions boundary_conditions;
    boundary_conditions.define(simParams().center, bparms,
                               amrex::DefaultGeometry(), nghost);

    amrex::Vector<amrex::BCRec> bcs(NUM_VARS);
    for (int icomp = 0; icomp < NUM_VARS; ++icomp)
    {
        auto &bc = bcs[icomp];
        for (amrex::OrientationIter oit; oit.isValid(); ++oit)
        {
            amrex::Orientation face = oit();
            const int idim          = face.coordDir();
            const int bctype = boundary_conditions.get_boundary_condition(face);
            if (amrex::DefaultGeometry().isPeriodic(idim))
            {
                bc.set(face, amrex::BCType::int_dir);
            }
            else if (bctype == BoundaryConditions::STATIC_BC ||
                     bctype == BoundaryConditions::SOMMERFELD_BC ||
                     bctype == BoundaryConditions::MIXED_BC)
            {
                bc.set(face, amrex::BCType::foextrap);
            }
            else if (bctype == BoundaryConditions::REFLECTIVE_BC)
            {
                int parity =
                    BoundaryConditions::get_state_var_parity(icomp, idim);
                if (parity == 1)
                {
                    bc.set(face, amrex::BCType::reflect_even);
                }
                else
                {
                    bc.set(face, amrex::BCType::reflect_odd);
                }
            }
            else if (bctype == BoundaryConditions::EXTRAPOLATING_BC)
            {
                amrex::Abort("xxxxx EXTRAPOLATING_BC todo");
            }
            else
            {
                amrex::Abort("Unknow BC type " + std::to_string(bctype));
            }
        }
    }

    amrex::StateDescriptor::BndryFunc bndryfunc(null_bc_fill);
    bndryfunc.setRunOnGPU(true); // Run the bc function on gpu.

    desc_lst.setComponent(zero_state_index, 0, StateVariables::names, bcs,
                          bndryfunc);
}

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
        4, zero_state_index, time, dt, iteration, ncycle,
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
        amrex::MultiFab &S_fine = fine_level.get_new_data(zero_state_index);
        amrex::MultiFab &S_crse = this->get_new_data(zero_state_index);
        amrex::Real t           = get_state_data(zero_state_index).curTime();

        amrex::IntVect ratio = parent->refRatio(lev);
        AMREX_ASSERT(ratio == 2 || ratio == 4);
        if (ratio == 2)
        {
            // Need to fill one ghost cell for the high-order interpolation
            // below
            FillPatch(fine_level, S_fine, 1, t, zero_state_index, 0,
                      S_fine.nComp());
        }

        FourthOrderInterpFromFineToCoarse(S_crse, 0, NUM_VARS, S_fine, ratio);
    }

    specificPostTimeStep();
}

void GRAMRLevel::post_regrid(int a_lbase, int a_new_finest)
{
    specific_post_regrid(a_lbase, a_new_finest);
}

void GRAMRLevel::post_init(amrex::Real /*stop_time*/)
{
    if (Level() == 0)
    {
        get_gramr_ptr()->set_restart_time(get_gramr_ptr()->cumTime());
    }
    specific_post_init();
}

void GRAMRLevel::post_restart()
{
    if (Level() == 0)
    {
        get_gramr_ptr()->set_restart_time(get_gramr_ptr()->cumTime());
    }
    specific_post_restart();
}

void GRAMRLevel::init(amrex::AmrLevel &old)
{
    BL_PROFILE("GRAMRLevel::init()");
    amrex::Real dt_new    = parent->dtLevel(level);
    amrex::Real cur_time  = old.get_state_data(zero_state_index).curTime();
    amrex::Real prev_time = old.get_state_data(zero_state_index).prevTime();
    amrex::Real dt_old    = cur_time - prev_time;
    setTimeLevel(cur_time, dt_old, dt_new);

    amrex::MultiFab &S_new = get_new_data(zero_state_index);
    FillPatch(old, S_new, 0, cur_time, zero_state_index, 0, S_new.nComp());
}

void GRAMRLevel::init()
{
    BL_PROFILE("GRAMRLevel::init()");
    amrex::Real dt = parent->dtLevel(level);
    const auto &coarse_state =
        parent->getLevel(level - 1).get_state_data(zero_state_index);
    amrex::Real cur_time  = coarse_state.curTime();
    amrex::Real prev_time = coarse_state.prevTime();
    amrex::Real dt_old =
        (cur_time - prev_time) /
        static_cast<amrex::Real>(parent->MaxRefRatio(level - 1));
    setTimeLevel(cur_time, dt_old, dt);

    amrex::MultiFab &S_new = get_new_data(zero_state_index);
    FillCoarsePatch(S_new, 0, cur_time, zero_state_index, 0, S_new.nComp());
}

void GRAMRLevel::errorEst(amrex::TagBoxArray &a_tag_box_array,
                          int /*a_clearval*/, int /*a_tagval*/,
                          amrex::Real /*a_time*/, int /*a_n_error_buf*/,
                          int /*a_ngrow*/)
{
    BL_PROFILE("GRAMRLevel::errorEst()");

    pre_tag_cells();

    // It is up to the derived class to use regrid_threshold in tag_cells()
    amrex::Real regrid_threshold = simParams().regrid_thresholds[Level()];
    tag_cells(a_tag_box_array, regrid_threshold);
}

void GRAMRLevel::writePlotFilePre(const std::string &a_dir, std::ostream &a_os)
{
    specific_pre_plotfile(a_dir, a_os);
}

void GRAMRLevel::writePlotFilePost(const std::string &a_dir, std::ostream &a_os)
{
    specific_post_plotfile(a_dir, a_os);
}

void GRAMRLevel::checkPointPre(const std::string &a_dir, std::ostream &a_os)
{
    specific_pre_checkpoint(a_dir, a_os);
}

void GRAMRLevel::checkPointPost(const std::string &a_dir, std::ostream &a_os)
{
    specific_post_checkpoint(a_dir, a_os);
}

bool GRAMRLevel::at_level_timestep_multiple(int a_level)
{
    // handle both the case a_level < Level() and a_level >= Level()
    int coarser_level     = std::min(a_level, Level());
    int finer_level       = std::max(a_level, Level());
    int finer_level_steps = get_gramr_ptr()->levelSteps(finer_level);

    // work out what the coarser level step number corresponds to on the finer
    // level
    int coarser_level_steps_at_finer_level =
        get_gramr_ptr()->levelSteps(coarser_level);

    for (int ilev = coarser_level + 1; ilev <= finer_level; ++ilev)
    {
        coarser_level_steps_at_finer_level *= get_gramr_ptr()->nCycle(ilev);
    }
    // finer_level_steps will be > coarser_level_steps
    return (finer_level_steps == coarser_level_steps_at_finer_level);
}