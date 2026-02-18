#include "WormholeLevel.hpp"
#include "CCZ4RHS.hpp"
#include "ChiTagger.hpp"
#include "Constraints.hpp"
#include "GRParmParse.hpp"
#include "PositiveChiAndLapse.hpp"
#include "SmallDataIO.hpp"
#include "TraceARemoval.hpp"
#include "Weyl4.hpp"
#include "WeylExtraction.hpp"
#include "WormholeInitialData.hpp"

#include <AMReX_Reduce.H>
#include <AMReX_Utility.H>

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
    BL_PROFILE("WormholeLevel::specificPostTimeStep");

    // --- Constraint norms (small-data output) ---
    // Compute volume-weighted L2 norms of Ham and Mom on level 0.
    if (simParams().calculate_constraint_norms && Level() == 0)
    {
        const amrex::Real time         = get_state_data(State_Type).curTime();
        const amrex::Real dt           = parent->dtLevel(0);
        const amrex::Real restart_time = get_gramr_ptr()->get_restart_time();
        const bool first_step          = (time == 0.0);

        // Fill ghosts for constraint calculation
        amrex::MultiFab &state_new = get_new_data(State_Type);
        FillPatch(*this, state_new, 2, time, State_Type, 0, state_new.nComp());

        // Compute constraints into a temporary MultiFab (Ham, Mom1, Mom2, Mom3)
        amrex::MultiFab cst(state_new.boxArray(), state_new.DistributionMap(), 4,
                            0);
        cst.setVal(0.0);
        Constraints::compute_mf(cst, 0, 4, state_new, Geom(), time, nullptr, 0);

        // Volume-weighted reductions on this level (constant cell volume on a level)
        const auto dx = Geom().CellSizeArray();
        const amrex::Real cell_vol = dx[0] * dx[1] * dx[2];

        amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpSum,
                         amrex::ReduceOpSum>
            reduce_ops;
        amrex::ReduceData<amrex::Real, amrex::Real, amrex::Real> reduce_data(
            reduce_ops);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        for (amrex::MFIter mfi(cst, amrex::TilingIfNotGPU()); mfi.isValid();
             ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            const auto arr       = cst.const_array(mfi);
            reduce_ops.eval(
                bx, reduce_data,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple
                {
                    const amrex::Real ham = arr(i, j, k, 0);
                    const amrex::Real m1  = arr(i, j, k, 1);
                    const amrex::Real m2  = arr(i, j, k, 2);
                    const amrex::Real m3  = arr(i, j, k, 3);
                    const amrex::Real mom2 = (m1 * m1 + m2 * m2 + m3 * m3);
                    return {ham * ham * cell_vol, mom2 * cell_vol, cell_vol};
                });
        }

        auto [sum_ham2, sum_mom2, sum_vol] = reduce_data.value();
        amrex::ParallelDescriptor::ReduceRealSum(sum_ham2);
        amrex::ParallelDescriptor::ReduceRealSum(sum_mom2);
        amrex::ParallelDescriptor::ReduceRealSum(sum_vol);

        const double L2_Ham =
            (sum_vol > 0.0) ? std::sqrt(sum_ham2 / sum_vol) : 0.0;
        const double L2_Mom =
            (sum_vol > 0.0) ? std::sqrt(sum_mom2 / sum_vol) : 0.0;

        // Build output directory: output_path + data_subpath
        GRParmParse pp;
        std::string output_path = "./";
        pp.load("output_path", output_path, std::string("./"));
        std::string data_subpath;
        pp.load("data_subpath", data_subpath, std::string(""));

        if (!output_path.empty() && output_path.back() != '/')
            output_path += "/";
        if (!data_subpath.empty() && data_subpath.back() != '/')
            data_subpath += "/";

        const std::string out_dir = output_path + data_subpath;
        if (!out_dir.empty())
        {
            amrex::UtilCreateDirectory(out_dir, 0755, false);
        }

        const std::string prefix = out_dir + "constraint_norms";

        SmallDataIO constraints_file(prefix, dt, time, restart_time,
                                     SmallDataIO::APPEND, first_step);
        constraints_file.remove_duplicate_time_data();
        if (first_step)
        {
            constraints_file.write_header_line({"L2_Ham", "L2_Mom"});
        }
        constraints_file.write_time_data_line({L2_Ham, L2_Mom});
    }
}