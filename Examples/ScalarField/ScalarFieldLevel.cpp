/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// General includes common to most GR problems
#include "ScalarFieldLevel.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "TraceARemoval.hpp"
// //#include "SixthOrderDerivatives.hpp"

// // For RHS update
#include "MatterCCZ4RHS.hpp"

// // For constraints calculation
// #include "NewMatterConstraints.hpp"

// // Problem specific includes
// //#include "ComputePack.hpp"
// //#include "GammaCalculator.hpp"

#include "InitialScalarData.hpp"
#include "Potential.hpp"
#include "ScalarField.hpp"
// //#include "SetValue.hpp"

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
    InitialScalarData gaussian(simParams().initial_params, dx[0]);

    amrex::MultiFab &state  = get_new_data(State_Type);
    auto const &state_array = state.arrays();

    amrex::ParallelFor(
        state, state.nGrowVect(),
        [=] AMREX_GPU_DEVICE(int box_ind, int i, int j, int k) noexcept
        { gaussian.compute(i, j, k, state_array[box_ind]); });
}

#ifdef AMREX_USE_HDF5
// Things to do before outputting a checkpoint file
void ScalarFieldLevel::prePlotLevel()
{
    // fillAllGhosts();
    // Potential potential(m_p.potential_params);
    // ScalarFieldWithPotential scalar_field(potential);
    // BoxLoops::loop(
    //     MatterConstraints<ScalarFieldWithPotential>(
    //         scalar_field, m_dx, m_p.G_Newton, c_Ham, Interval(c_Mom, c_Mom)),
    //     m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
}
#endif

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
    Potential potential(simParams().potential_params);
    ScalarFieldWithPotential scalar_field(potential);

    // Calculate CCZ4 right hand side
    if (simParams().max_spatial_derivative_order == 4)
    {
        MatterCCZ4RHS<ScalarFieldWithPotential, MovingPunctureGauge,
                      FourthOrderDerivatives>
            matter_ccz4_rhs(scalar_field, simParams().ccz4_params,
                            Geom().CellSize(0), simParams().sigma,
                            simParams().formulation, simParams().G_Newton);
        amrex::ParallelFor(
            a_rhs,
            [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) {
                matter_ccz4_rhs.compute(i, j, k, rhs_arrs[box_no],
                                        soln_c_arrs[box_no]);
            });
    }
    else if (simParams().max_spatial_derivative_order == 6)
    {
        amrex::Abort("xxxxx max_spatial_derivative_order == 6 todo");
#if 0
        MatterCCZ4RHS<ScalarFieldWithPotential, MovingPunctureGauge, SixthOrderDerivatives>
	  matter_ccz4_rhs(scalar_field, simParams().ccz4_params, Geom().CellSize(0), simParams().sigma,
			  simParams().formulation, simParams().G_Newton);
        amrex::ParallelFor(a_rhs,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k)
        {
            amrex::CellData<amrex::Real const> state = soln_c_arrs[box_no].cellData(i,j,k);
            amrex::CellData<amrex::Real> rhs = rhs_arrs[box_no].cellData(i,j,k);
            matter_ccz4_rhs.compute(i,j,k,rhs_arrs[box_no], soln_c_arrs[box_no]);
        });
#endif
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

void ScalarFieldLevel::preTagCells()
{
    // we don't need any ghosts filled for the fixed grids tagging criterion
    // used here so don't fill any
}

void ScalarFieldLevel::errorEst(amrex::TagBoxArray &tagging_criterion,
                                int clearval, int tagval, amrex::Real time,
                                int n_error_buf, int ngrow)

{
    BL_PROFILE("ScalarFieldLevel::errorEst()");
    // BoxLoops::loop(
    //     FixedGridsTaggingCriterion(m_dx, m_level, 2.0 * m_p.L, m_p.center),
    //     current_state, tagging_criterion);
}
