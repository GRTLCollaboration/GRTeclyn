/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "GRAMRLevel.hpp"

struct GRAMRBCFill
{
    AMREX_GPU_DEVICE
    void operator() (const amrex::IntVect& /*iv*/,
                     amrex::Array4<amrex::Real> const& /*dest*/,
                     const int /*dcomp*/, const int /*numcomp*/,
                     amrex::GeometryData const& /*geom*/,
                     const amrex::Real /*time*/,
                     const amrex::BCRec* /*bcr*/, const int /*bcomp*/,
                     const int /*orig_comp*/) const
    {
        // xxxxx BCFill todo
    }
};

void gramr_bc_fill(amrex::Box const& bx, amrex::FArrayBox& data,
                   const int dcomp, const int numcomp,
                   amrex::Geometry const& geom, const amrex::Real time,
                   const amrex::Vector<amrex::BCRec>& bcr, const int bcomp,
                   const int scomp)
{
    amrex::GpuBndryFuncFab<GRAMRBCFill> bndry_func(GRAMRBCFill{});
    bndry_func(bx,data,dcomp,numcomp,geom,time,bcr,bcomp,scomp);
}

void GRAMRLevel::variableSetUp()
{
    desc_lst.addDescriptor(State_Type,amrex::IndexType::TheCellType(),
                           amrex::StateDescriptor::Point,
                           simParams().num_ghosts,
                           NUM_VARS,
                           &amrex::cell_quartic_interp);

    int lo_bc[BL_SPACEDIM];
    int hi_bc[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; ++i) {
        if (amrex::DefaultGeometry().isPeriodic(i)) {
            lo_bc[i] = hi_bc[i] = amrex::BCType::int_dir;
        } else {
            lo_bc[i] = hi_bc[i] = amrex::BCType::ext_dir;
        }
    }

    amrex::Vector<amrex::BCRec> bcs(NUM_VARS);
    amrex::Vector<std::string> name(NUM_VARS);
    for (int i = 0; i < NUM_VARS; ++i) {
        bcs[i] = amrex::BCRec(lo_bc, hi_bc);
        name[i] = UserVariables::variable_names[i];
    }

    amrex::StateDescriptor::BndryFunc bndryfunc(gramr_bc_fill);
    bndryfunc.setRunOnGPU(true);  // Run the bc function on gpu.

    desc_lst.setComponent(State_Type, 0, name, bcs, bndryfunc);
}

void GRAMRLevel::variableCleanUp()
{
    desc_lst.clear();
}

GRAMRLevel::GRAMRLevel() {}

GRAMRLevel::GRAMRLevel(amrex::Amr& papa, int lev,
                       const amrex::Geometry& a_geom,
                       const amrex::BoxArray& ba,
                       const amrex::DistributionMapping& dm,
                       amrex::Real time)
    : amrex::AmrLevel(papa,lev,a_geom,ba,dm,time)
{
    m_num_ghosts = simParams().num_ghosts;
    m_boundaries.define(simParams().center, simParams().boundary_params,
                        a_geom, m_num_ghosts);
}

GRAMRLevel::~GRAMRLevel() {}

SimulationParameters const& GRAMRLevel::simParams ()
{
    return GRAMR::get_simulation_parameters();
}

void GRAMRLevel::computeInitialDt (int finest_level, int /*sub_cycle*/,
                                   amrex::Vector<int>& /*n_cycle*/,
                                   const amrex::Vector<amrex::IntVect>&/*ref_ratio*/,
                                   amrex::Vector<amrex::Real>& dt_level,
                                   amrex::Real /*stop_time*/)
{
    // Level 0 will do it for all levels
    if (Level() == 0) {
        double dt_multiplier = simParams().dt_multiplier;
        for (int i = 0; i <= finest_level; ++i) {
            dt_level[i] = dt_multiplier * parent->Geom(i).CellSize(0);
        }
    }
}

void GRAMRLevel::computeNewDt (int finest_level, int /*sub_cycle*/,
                               amrex::Vector<int>& /*n_cycle*/,
                               const amrex::Vector<amrex::IntVect>& /*ref_ratio*/,
                               amrex::Vector<amrex::Real>& dt_min,
                               amrex::Vector<amrex::Real>& dt_level,
                               amrex::Real /*stop_time*/, int /*post_regrid_flag*/)
{
    // This is called at the end of a coarse time step
    // Level 0 will do it for all levels
    if (Level() == 0) {
        double dt_multiplier = simParams().dt_multiplier;
        for (int i = 0; i <= finest_level; ++i) {
            dt_min[i] = dt_level[i] = dt_multiplier * parent->Geom(i).CellSize(0);
        }
    }
}

amrex::Real GRAMRLevel::advance (amrex::Real time, amrex::Real dt,
                                 int iteration, int ncycle)
{
    for (int k = 0; k < NUM_STATE_TYPE; k++) {
        state[k].allocOldData();
        state[k].swapTimeLevels(dt);
    }

    amrex::AmrLevel::RK(4, State_Type, time, dt, iteration, ncycle,
        [&] (int /*stage*/, amrex::MultiFab& dSdt,
             amrex::MultiFab const& S,
             amrex::Real t, amrex::Real /*dtsub*/)
        {
            specificEvalRHS(const_cast<amrex::MultiFab&>(S), dSdt, t);
        },
        [&] (int /*stage*/, amrex::MultiFab& S) {
            specificUpdateODE(S);
        });

    specificAdvance();

    return dt;
}

void GRAMRLevel::post_timestep (int /*iteration*/)
{
    const int lev = Level();
    if (lev < parent->finestLevel()) {
        auto& fine_level = parent->getLevel(Level()+1);
        amrex::MultiFab & S_fine = fine_level.get_new_data(State_Type);
        amrex::MultiFab & S_crse =      this->get_new_data(State_Type);
        amrex::Real t = get_state_data(State_Type).curTime();

        amrex::IntVect ratio = parent->refRatio(lev);
        AMREX_ASSERT(ratio == 2 || ratio == 4);
        if (ratio == 2) {
            // Need to fill one ghost cell for the high-order interpolation below
            FillPatch(fine_level, S_fine, 1, t, State_Type, 0, S_fine.nComp());
        }

        FourthOrderInterpFromFineToCoarse(S_crse, 0, 2, S_fine, ratio);
    }

    specificPostTimeStep();
}

void GRAMRLevel::post_regrid (int lbase, int new_finest)
{
    amrex::ignore_unused(lbase, new_finest);
}

void GRAMRLevel::post_init (amrex::Real /*stop_time*/)
{
    // Don't we need to do anything here
}

void GRAMRLevel::init (amrex::AmrLevel &old)
{
    amrex::ignore_unused(old);
    amrex::Abort("xxxxx GRAMRLevel::init(old) todo");
}

void GRAMRLevel::init ()
{
    amrex::Abort("xxxxx GRAMRLevel::init() todo");
}

void GRAMRLevel::errorEst (amrex::TagBoxArray& tb, int clearval, int tagval,
                           amrex::Real time, int n_error_buf, int ngrow)
{
    amrex::ignore_unused(tb,clearval,tagval,time,n_error_buf,ngrow);
    amrex::Abort("xxxxx GRAMRLevel::errorEst todo");
}
