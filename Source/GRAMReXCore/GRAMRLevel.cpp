/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "GRAMRLevel.hpp"

amrex::Vector<std::string> GRAMRLevel::plot_constraints;

struct GRAMRBCFill
{
    AMREX_GPU_DEVICE
    void operator()(const amrex::IntVect & /*iv*/,
                    amrex::Array4<amrex::Real> const & /*dest*/,
                    const int /*dcomp*/, const int /*numcomp*/,
                    amrex::GeometryData const & /*geom*/,
                    const amrex::Real /*time*/, const amrex::BCRec * /*bcr*/,
                    const int /*bcomp*/, const int /*orig_comp*/) const
    {
        // We don't have any external Dirichlet BC
    }
};

void gramr_bc_fill(amrex::Box const &bx, amrex::FArrayBox &data,
                   const int dcomp, const int numcomp,
                   amrex::Geometry const &geom, const amrex::Real time,
                   const amrex::Vector<amrex::BCRec> &bcr, const int bcomp,
                   const int scomp)
{
    amrex::GpuBndryFuncFab<GRAMRBCFill> bndry_func(GRAMRBCFill{});
    bndry_func(bx, data, dcomp, numcomp, geom, time, bcr, bcomp, scomp);
}

void GRAMRLevel::variableSetUp()
{
    const int nghost = simParams().num_ghosts;
    desc_lst.addDescriptor(State_Type, amrex::IndexType::TheCellType(),
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
            else if (bctype == BoundaryConditions::STATIC_BC)
            {
                bc.set(face, amrex::BCType::foextrap);
            }
            else if (bctype == BoundaryConditions::SOMMERFELD_BC)
            {
                bc.set(face, amrex::BCType::foextrap);
            }
            else if (bctype == BoundaryConditions::REFLECTIVE_BC)
            {
                int parity = boundary_conditions.get_var_parity(
                    icomp, idim, VariableType::evolution);
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
            else if (bctype == BoundaryConditions::MIXED_BC)
            {
                bc.set(face, amrex::BCType::foextrap); // xxxxx todo mixed BC
            }
            else
            {
                amrex::Abort("Unknow BC type " + std::to_string(bctype));
            }
        }
    }

    amrex::Vector<std::string> name(NUM_VARS);
    for (int i = 0; i < NUM_VARS; ++i)
    {
        name[i] = UserVariables::variable_names[i];
    }

    amrex::StateDescriptor::BndryFunc bndryfunc(gramr_bc_fill);
    bndryfunc.setRunOnGPU(true); // Run the bc function on gpu.

    desc_lst.setComponent(State_Type, 0, name, bcs, bndryfunc);

    amrex::ParmParse pp("amr");
    if (pp.contains("derive_plot_vars"))
    {
        std::vector<std::string> names;
        pp.getarr("derive_plot_vars", names);

        // Constraints
        auto r =
            std::find_if(names.begin(), names.end(),
                         [](std::string const &s)
                         { return std::string("ham") == amrex::toLower(s); });
        if (r != names.end())
        {
            names.erase(r);
            plot_constraints.push_back("Ham");
        }
        //
        r = std::find_if(names.begin(), names.end(),
                         [](std::string const &s)
                         { return std::string("mom") == amrex::toLower(s); });
        if (r != names.end())
        {
            names.erase(r);
            plot_constraints.push_back("Mom1");
            plot_constraints.push_back("Mom2");
            plot_constraints.push_back("Mom3");
        }
        //
        if (!plot_constraints.empty())
        {
            names.push_back("constraints");
            pp.addarr("derive_plot_vars", names);

            derive_lst.add(
                "constraints", amrex::IndexType::TheCellType(),
                plot_constraints.size(), plot_constraints,
                amrex::DeriveFuncFab(), // null function because we won't use
                                        // it.
                [=](amrex::Box const &b) { return amrex::grow(b, nghost); },
                &amrex::cell_quartic_interp);
            derive_lst.addComponent("constraints", desc_lst, State_Type, 0,
                                    NUM_VARS);
        }
    }
}

void GRAMRLevel::variableCleanUp()
{
    desc_lst.clear();
    derive_lst.clear();
}

GRAMRLevel::GRAMRLevel() {}

GRAMRLevel::GRAMRLevel(amrex::Amr &papa, int lev, const amrex::Geometry &a_geom,
                       const amrex::BoxArray &ba,
                       const amrex::DistributionMapping &dm, amrex::Real time)
    : amrex::AmrLevel(papa, lev, a_geom, ba, dm, time)
{
    m_num_ghosts = simParams().num_ghosts;
    m_boundaries.define(simParams().center, simParams().boundary_params, a_geom,
                        m_num_ghosts);
}

GRAMRLevel::~GRAMRLevel() {}

SimulationParameters const &GRAMRLevel::simParams()
{
    return GRAMR::get_simulation_parameters();
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
    for (int k = 0; k < NUM_STATE_TYPE; k++)
    {
        state[k].allocOldData();
        state[k].swapTimeLevels(dt);
    }

    amrex::AmrLevel::RK(
        4, State_Type, time, dt, iteration, ncycle,
        [&](int /*stage*/, amrex::MultiFab &dSdt, amrex::MultiFab const &S,
            amrex::Real t, amrex::Real /*dtsub*/)
        {
            specificEvalRHS(const_cast<amrex::MultiFab &>(S), dSdt, t);
            m_boundaries.apply_sommerfeld_boundaries(dSdt, S);
        },
        [&](int /*stage*/, amrex::MultiFab &S) { specificUpdateODE(S); });

    specificAdvance();

    return dt;
}

void GRAMRLevel::post_timestep(int /*iteration*/)
{
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
    // xxxxx Do we need to do anything after the initializaion?
}

void GRAMRLevel::init(amrex::AmrLevel &old)
{
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
    amrex::Real dt = parent->dtLevel(level);
    auto const &coarse_state =
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
    if (!plot_constraints.empty())
    {
        auto &S = get_new_data(State_Type);
        FillPatch(*this, S, S.nGrow(), get_state_data(State_Type).curTime(),
                  State_Type, 0, S.nComp()); // xxxxx Do we need all components?
    }
}

void GRAMRLevel::writePlotFilePost(const std::string & /*dir*/,
                                   std::ostream & /*os*/)
{
    m_is_writing_plotfile = false;
}
