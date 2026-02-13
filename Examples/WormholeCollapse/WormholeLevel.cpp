#include "WormholeLevel.hpp"
#include "CCZ4RHS.hpp"
#include "ChiTagger.hpp"
#include "Constraints.hpp"
#include "PositiveChiAndLapse.hpp"
#include "TraceARemoval.hpp"
#include "Weyl4.hpp"
#include "WeylExtraction.hpp"
#include "WormholeInitialData.hpp"

void WormholeLevel::variableSetUp()
{
    BL_PROFILE("WormholeLevel::variableSetUp()");
    stateVariableSetUp();
    Constraints::set_up(State_Type);
    Weyl4::set_up(State_Type);
}

void WormholeLevel::specificAdvance()
{
    amrex::MultiFab &S_new = get_new_data(State_Type);
    const auto &arrs       = S_new.arrays();

    // Enforce algebraic constraints
    amrex::ParallelFor(S_new,
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       {
                           amrex::CellData<amrex::Real> cell =
                               arrs[box_no].cellData(i, j, k);
                           TraceARemoval()(cell);
                           PositiveChiAndLapse()(cell);
                       });
}

void WormholeLevel::initData()
{
    BL_PROFILE("WormholeLevel::initData");

    // Set up the compute class for the Wormhole initial data
    WormholeInitialData wormhole(simParams().wormhole_params,
                                 Geom().CellSize(0));

    amrex::MultiFab &state = get_new_data(State_Type);
    const auto &arrs       = state.arrays();

    amrex::ParallelFor(state, state.nGrowVect(),
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       {
                           // Initialize all to zero first
                           amrex::CellData<amrex::Real> cell =
                               arrs[box_no].cellData(i, j, k);
                           for (int n = 0; n < cell.nComp(); ++n)
                           {
                               cell[n] = 0.;
                           }
                           // Calculate Wormhole metrics
                           wormhole.compute(i, j, k, arrs[box_no]);
                       });

    amrex::Gpu::streamSynchronize();
}

void WormholeLevel::specificEvalRHS(amrex::MultiFab &a_soln,
                                    amrex::MultiFab &a_rhs,
                                    const double /*a_time*/)
{
    BL_PROFILE("WormholeLevel::specificEvalRHS()");
    const auto &soln_arrs   = a_soln.arrays();
    const auto &soln_c_arrs = a_soln.const_arrays();
    const auto &rhs_arrs    = a_rhs.arrays();

    // Enforce algebraic constraints pre-RHS
    amrex::ParallelFor(a_soln, a_soln.nGrowVect(),
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       {
                           amrex::CellData<amrex::Real> cell =
                               soln_arrs[box_no].cellData(i, j, k);
                           TraceARemoval()(cell);
                           PositiveChiAndLapse()(cell);
                       });

    // Calculate CCZ4 Right Hand Side (Vacuum Einstein Equations)
    CCZ4RHS<MovingPunctureGauge, FourthOrderDerivatives> ccz4rhs(
        simParams().ccz4_params, Geom().CellSize(0), simParams().sigma,
        simParams().formulation);

    amrex::ParallelFor(
        a_rhs, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
        { ccz4rhs.compute(i, j, k, rhs_arrs[box_no], soln_c_arrs[box_no]); });

    amrex::Gpu::streamSynchronize();
}

void WormholeLevel::specificUpdateODE(amrex::MultiFab &a_soln)
{
    const auto &soln_arrs = a_soln.arrays();
    amrex::ParallelFor(a_soln, amrex::IntVect(0),
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       {
                           amrex::CellData<amrex::Real> cell =
                               soln_arrs[box_no].cellData(i, j, k);
                           TraceARemoval()(cell);
                       });

    amrex::Gpu::streamSynchronize();
}

void WormholeLevel::pre_tag_cells()
{
    amrex::MultiFab &state_new = get_new_data(State_Type);
    const auto cur_time        = get_state_data(State_Type).curTime();
    FillPatch(*this, state_new, 2, cur_time, State_Type, c_chi, 1);
}

void WormholeLevel::tag_cells(amrex::TagBoxArray &a_tag_box_array,
                              amrex::Real a_regrid_threshold)
{
    BL_PROFILE("WormholeLevel::tag_cells()");
    amrex::MultiFab &state_new = get_new_data(State_Type);
    const auto &tag_arrs       = a_tag_box_array.arrays();
    const auto &state_new_arrs = state_new.const_arrays();

    // Tag based on gradients of Chi (Conformal factor)
    // This will naturally tag the throat region where curvature is high
    ChiTagger chi_tagger(Geom().CellSize(0), a_regrid_threshold);

    amrex::ParallelFor(state_new, amrex::IntVect(0),
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       {
                           const auto &tags_arr  = tag_arrs[box_no];
                           const auto &state_arr = state_new_arrs[box_no];
                           chi_tagger(i, j, k, tags_arr, state_arr);
                       });
    amrex::Gpu::streamSynchronize();
}

void WormholeLevel::specificPostTimeStep()
{
    // Implement diagnostics here (e.g., Weyl4 extraction)
    // This is similar to BinaryBH but stripped of PunctureTracker
}