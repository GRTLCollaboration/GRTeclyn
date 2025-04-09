/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// General includes common to most GR problems
#include "ScalarFieldLevel.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "TraceARemoval.hpp"

// For RHS update
#include "CCZ4RHSWithMatter.hpp"

// For constraints calculation
#include "Constraints.hpp"
#include "ConstraintsWithMatter.hpp"
#include "Weyl4WithMatter.hpp"

// For tagging cells
#include "ChiExtractionTagger.hpp"
#include "FixedGridsTagger.hpp"

// Problem specific includes
#include "InitialScalarData.hpp"
#include "Potential.hpp"
#include "ScalarField.hpp"

void ScalarFieldLevel::variableSetUp()
{
    BL_PROFILE("ScalarFieldLevel::variableSetUp()");

    // Set up the state variables
    stateVariableSetUp();

    const int nghost = simParams().num_ghosts;

    // // Add the constraints to the derive list
    derive_lst.add(
        "constraints", amrex::IndexType::TheCellType(),
        static_cast<int>(Constraints::var_names.size()), Constraints::var_names,
        amrex::DeriveFuncFab(), // null function because we won't use
                                // it.
        [=](const amrex::Box &box) { return amrex::grow(box, nghost); },
        &amrex::cell_quartic_interp);

    derive_lst.addComponent("constraints", desc_lst, State_Type, 0, NUM_VARS);

    // Add Weyl4 to the derive list
    derive_lst.add(
        "Weyl4", amrex::IndexType::TheCellType(),
        static_cast<int>(Weyl4::var_names.size()), Weyl4::var_names,
        amrex::DeriveFuncFab(), // null function because we won't use it
        [=](const amrex::Box &box) { return amrex::grow(box, nghost); },
        &amrex::cell_quartic_interp);

    derive_lst.addComponent("Weyl4", desc_lst, State_Type, 0, NUM_VARS);
}

// Things to do at each advance step, after the RK4 is calculated
void ScalarFieldLevel::specificAdvance()
{
    BL_PROFILE("ScalarFieldLevel::specificAdvance");
    // Enforce trace free A_ij and positive chi and alpha
    amrex::MultiFab &S_new = get_new_data(State_Type);
    const auto &arrs       = S_new.arrays();

    // Enforce the trace free A_ij condition and positive chi and alpha
    amrex::ParallelFor(S_new,
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       {
                           amrex::CellData<amrex::Real> cell =
                               arrs[box_no].cellData(i, j, k);
                           TraceARemoval()(cell);
                           PositiveChiAndAlpha()(cell);
                       });

    // Check for nan's
    if (simParams().nan_check)
    {
        if (S_new.contains_nan(0, S_new.nComp(), amrex::IntVect(0), true))
        {
            amrex::Abort("NaN in specificAdvance");
        }
    }
}

// Initial data for field and metric variables
void ScalarFieldLevel::initData()
{
    BL_PROFILE("ScalarFieldLevel::initData");
    if (m_verbosity)
        amrex::Print() << "ScalarFieldLevel::initialData " << Level()
                       << std::endl;

    const auto dx = geom.CellSizeArray();
    InitialScalarData gaussian_pulse(simParams().initial_params, dx[0]);

    amrex::MultiFab &state  = get_new_data(State_Type);
    auto const &state_array = state.arrays();

    amrex::ParallelFor(
        state, state.nGrowVect(),
        [=] AMREX_GPU_DEVICE(int box_ind, int i, int j, int k) noexcept
        {
            amrex::CellData<amrex::Real> cell =
                state_array[box_ind].cellData(i, j, k);
            for (int n = 0; n < cell.nComp(); ++n)
            {
                cell[n] = 0.;
            }

            gaussian_pulse.compute(i, j, k, state_array[box_ind]);
        });

    if (simParams().nan_check)
    {
        if (state.contains_nan(0, state.nComp(), amrex::IntVect(0), true))
        {
            amrex::Abort("NaN in initData");
        }
    }
}

// Things to do in RHS update, at each RK4 step
void ScalarFieldLevel::specificEvalRHS(amrex::MultiFab &a_soln,
                                       amrex::MultiFab &a_rhs,
                                       const double a_time)
{
    BL_PROFILE("ScalarFieldLevel::specificEvalRHS()");

    const auto &soln_arrs   = a_soln.arrays();
    const auto &soln_c_arrs = a_soln.const_arrays();
    const auto &rhs_arrs    = a_rhs.arrays();

    // Enforce positive chi and alpha and trace free A
    amrex::ParallelFor(a_soln, a_soln.nGrowVect(),
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       {
                           amrex::CellData<amrex::Real> cell =
                               soln_arrs[box_no].cellData(i, j, k);
                           TraceARemoval()(cell);
                           PositiveChiAndAlpha()(cell);
                       });

    // Calculate MatterCCZ4 right hand side with matter_t = ScalarField
    ScalarFieldWithPotential scalar_field;

    // Calculate CCZ4 right hand side
    if (simParams().max_spatial_derivative_order == 4)
    {
        CCZ4RHSWithMatter<ScalarFieldWithPotential,
                          MovingPunctureGaugeWithMatter, FourthOrderDerivatives>
            ccz4_rhs_with_matter(simParams().ccz4_params, Geom().CellSize(0),
                                 simParams().sigma, simParams().formulation,
                                 simParams().G_Newton);
        amrex::ParallelFor(
            a_rhs,
            // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
            [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
            {
                ccz4_rhs_with_matter.compute(i, j, k, rhs_arrs[box_no],
                                             soln_c_arrs[box_no]);
            });
    }
    else if (simParams().max_spatial_derivative_order == 6)
    {
        amrex::Abort("xxxxx max_spatial_derivative_order == 6 todo");
#if 0
        CCZ4RHSWithMatter<ScalarFieldWithPotential, MovingPunctureGaugeWithMatter, SixthOrderDerivatives>
	  ccz4_rhs_with_matter(simParams().ccz4_params, Geom().CellSize(0), simParams().sigma,
			  simParams().formulation, simParams().G_Newton);
        amrex::ParallelFor(a_rhs,
        //NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k)
        {
            amrex::CellData<amrex::Real const> state = soln_c_arrs[box_no].cellData(i,j,k);
            amrex::CellData<amrex::Real> rhs = rhs_arrs[box_no].cellData(i,j,k);
            ccz4_rhs_with_matter.compute(i,j,k,rhs_arrs[box_no], soln_c_arrs[box_no]);
        });
#endif
    }

    if (simParams().nan_check)
    {
        if (a_soln.contains_nan(0, a_soln.nComp(), amrex::IntVect(0), true))
        {

            for (int i = 0; i < NUM_VARS; i++)
            {
                if (a_soln.contains_nan(i, 1, amrex::IntVect(0), true))
                    amrex::Print() << "Nan found in "
                                   << StateVariables::names[i] << std::endl;
            }

            amrex::Abort("NaN in specificUpdateRHS");
        }
    }
}

// Things to do at ODE update, after soln + rhs
void ScalarFieldLevel::specificUpdateODE(amrex::MultiFab &a_soln)
{
    BL_PROFILE("ScalarFieldLevel::specificUpdateODE()");
    // Enforce the trace free A_ij condition
    const auto &soln_arrs = a_soln.arrays();
    amrex::ParallelFor(a_soln, amrex::IntVect(0), // zero ghost cells
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       {
                           amrex::CellData<amrex::Real> cell =
                               soln_arrs[box_no].cellData(i, j, k);
                           TraceARemoval()(cell);
                       });
}

void ScalarFieldLevel::pre_tag_cells()
{
    // we don't need any ghosts filled for the fixed grids tagging criterion
    // used here so don't fill any
    amrex::MultiFab &state_new = get_new_data(State_Type);
    const auto curr_time       = get_state_data(State_Type).curTime();

    const int nghost = 2;
    const int ncomp  = 1;

    FillPatch(*this, state_new, nghost, curr_time, State_Type, c_chi, ncomp);
}

void ScalarFieldLevel::tag_cells(amrex::TagBoxArray &a_tag_box_array,
                                 amrex::Real a_regrid_threshold)

{
    BL_PROFILE("ScalarFieldLevel::tag_cells()");

    amrex::MultiFab &state_new = get_new_data(State_Type);
    const auto curr_time       = get_state_data(State_Type).curTime();

    const auto &simpar = simParams();

    const amrex::Real sim_dx = Geom().CellSize(0);

    // Check for reflective BCs but there's only ever one short/reflected side
    const amrex::Real max_length =
        std::max(Geom().ProbLength(0), Geom().ProbLength(1));

    std::array<double, AMREX_SPACEDIM> sim_center =
        simpar.initial_params.center;

    const auto &tag_arrs       = a_tag_box_array.arrays();
    const auto &state_new_arrs = state_new.const_arrays();

    // FixedGridsTagger tagger(sim_dx, Level(), max_length, sim_center);
    amrex::Real threshold = simpar.regrid_thresholds[Level()];
    ChiExtractionTagger tagger(Geom().CellSize(0), Level(), threshold,
                               simParams().extraction_params,
                               simParams().activate_extraction);

    amrex::ParallelFor(state_new, amrex::IntVect(0),
                       // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       {
                           tagger(i, j, k, tag_arrs[box_no],
                                  state_new_arrs[box_no]);
                           // tagger(i, j, k, tag_arrs[box_no]);
                       });

    amrex::Gpu::streamSynchronize();
}

void ScalarFieldLevel::derive(const std::string &name, amrex::Real time,
                              amrex::MultiFab &multifab, int dcomp)
{
    BL_PROFILE("ScalarFieldLevel::derive()");

    BL_ASSERT(dcomp < multifab.nComp());

    const int num_ghosts = multifab.nGrow();

    const amrex::DeriveRec *rec = derive_lst.get(name);
    if (rec != nullptr)
    {
        int state_idx, derive_scomp, derive_ncomp;

        // we only have one state so state_idx will be State_Type = 0
        rec->getRange(0, state_idx, derive_scomp, derive_ncomp);

        // work out how many extra ghost cells we need
        const amrex::BoxArray &src_ba = state[state_idx].boxArray();

        int num_extra_ghosts = num_ghosts;
        {
            amrex::Box box0   = src_ba[0];
            amrex::Box box1   = rec->boxMap()(box0);
            num_extra_ghosts += box0.smallEnd(0) - box1.smallEnd(0);
        }

        // Make a Multifab with enough extra ghosts to calculated derived
        // quantity. For now use NUM_VARS in case the enum mapping loads more
        // vars than is actually needed
        amrex::MultiFab src_mf(src_ba, dmap, NUM_VARS, num_extra_ghosts,
                               amrex::MFInfo(), *m_factory);

        // Fill the multifab with the needed state data including the ghost
        // cells
        FillPatch(*this, src_mf, num_extra_ghosts, time, state_idx,
                  derive_scomp, derive_ncomp);

        const auto &src_arrays = src_mf.const_arrays();

        ScalarFieldWithPotential scalar_field;

        if (name == "constraints")
        {
            const auto &out_arrays = multifab.arrays();
            int iham               = dcomp;
            Interval imom = Interval(dcomp + 1, dcomp + AMREX_SPACEDIM);
            ConstraintsWithMatter<ScalarFieldWithPotential> constraints(
                Geom().CellSize(0), simParams().G_Newton, iham, imom);
            amrex::ParallelFor(
                multifab, multifab.nGrowVect(),
                [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
                {
                    constraints.compute(i, j, k, out_arrays[box_no],
                                        src_arrays[box_no]);
                });
        }
        else if (name == "Weyl4")
        {
            const auto &out_arrays = multifab.arrays();

            Weyl4WithMatter<ScalarFieldWithPotential> weyl4(
                simParams().extraction_params.center, Geom().CellSize(0), dcomp,
                simParams().formulation, simParams().G_Newton);

            amrex::ParallelFor(
                multifab, multifab.nGrowVect(),
                [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
                {
                    weyl4.compute(i, j, k, out_arrays[box_no],
                                  src_arrays[box_no]);
                });
        }
        else
        {
            amrex::Abort("Unknown derived variable");
        }
    }
    else
    {
        amrex::Abort("Unknown derived variable");
    }
    amrex::Gpu::streamSynchronize();
}
