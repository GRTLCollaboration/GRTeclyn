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
    Constraints::set_up(state_index);
    Weyl4::set_up(state_index);
}

void WormholeLevel::specificAdvance()
{
    amrex::MultiFab &S_new = get_new_data(state_index);
    const auto &arrs       = S_new.arrays();
    TraceARemoval trace_A_removal;
    PositiveChiAndLapse positive_chi_lapse;

    // Enforce algebraic constraints
    amrex::ParallelFor(S_new,
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       {
                           trace_A_removal(i, j, k, arrs[box_no]);
                           positive_chi_lapse(i, j, k, arrs[box_no]);
                       });
}

void WormholeLevel::initData()
{
    BL_PROFILE("WormholeLevel::initData");

    // Set up the compute class for the Wormhole initial data
    WormholeInitialData wormhole(simParams().wormhole_params,
                                 Geom().CellSize(0));

    amrex::MultiFab &state = get_new_data(state_index);
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
                                    const double a_time)
{
    BL_PROFILE("WormholeLevel::specificEvalRHS()");
    const auto &soln_arrs   = a_soln.arrays();
    const auto &soln_c_arrs = a_soln.const_arrays();
    const auto &rhs_arrs    = a_rhs.arrays();
    TraceARemoval trace_A_removal;
    PositiveChiAndLapse positive_chi_lapse;

    // Enforce algebraic constraints pre-RHS
    amrex::ParallelFor(a_soln, a_soln.nGrowVect(),
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       {
                           trace_A_removal(i, j, k, soln_arrs[box_no]);
                           positive_chi_lapse(i, j, k, soln_arrs[box_no]);
                       });

    // Calculate CCZ4 Right Hand Side (Vacuum Einstein Equations)
    CCZ4RHS<MovingPunctureGauge, FourthOrderDerivatives> ccz4rhs(
        simParams().ccz4_params, Geom().CellSize(0), simParams().sigma,
        simParams().formulation);

    amrex::ParallelFor(
        a_rhs, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
        { ccz4rhs(i, j, k, rhs_arrs[box_no], soln_c_arrs[box_no]); });

    const auto delayed_start =
        simParams().wormhole_params.delayed_kick_start;
    const auto delayed_duration =
        simParams().wormhole_params.delayed_kick_duration;
    const auto delayed_monopole =
        simParams().wormhole_params.delayed_kick_monopole_amplitude;
    const auto delayed_quadrupole =
        simParams().wormhole_params.delayed_kick_quadrupole_amplitude;
    const auto delayed_width =
        simParams().wormhole_params.delayed_kick_width;

    if (delayed_start >= 0.0 && delayed_duration > 0.0 &&
        ((delayed_monopole != 0.0) || (delayed_quadrupole != 0.0)) &&
        delayed_width > 0.0 && a_time >= delayed_start &&
        a_time <= delayed_start + delayed_duration)
    {
        const double x = (a_time - delayed_start) / delayed_duration;
        // Smooth pulse with zero endpoints. The time integral over the window is
        // unity, so the delayed amplitudes correspond approximately to the
        // injected change in K.
        const double temporal_weight =
            (1.0 - std::cos(2.0 * M_PI * x)) / delayed_duration;
        const double sig2 = delayed_width * delayed_width;
        const auto center = simParams().wormhole_params.grid_center;
        const auto centerA = simParams().wormhole_params.centerA;
        const auto centerB = simParams().wormhole_params.centerB;
        const int metric_type = simParams().wormhole_params.metric_type;
        const double dx = Geom().CellSize(0);

        amrex::ParallelFor(
            a_rhs, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
            {
                Coordinates coords(amrex::IntVect(AMREX_D_DECL(i, j, k)), dx,
                                   center);
                const amrex::Real xloc = coords.x;
                const amrex::Real yloc = coords.y;
                const amrex::Real zloc = coords.z;

                const amrex::Real dxA = xloc - centerA[0];
                const amrex::Real dyA = yloc - centerA[1];
                const amrex::Real dzA = zloc - centerA[2];
                const amrex::Real rA2 = dxA * dxA + dyA * dyA + dzA * dzA;
                const amrex::Real rA2_safe = amrex::max(rA2, amrex::Real(1.0e-24));
                const amrex::Real y20A = 3.0 * dzA * dzA / rA2_safe - 1.0;
                const amrex::Real kickA =
                    (delayed_monopole + delayed_quadrupole * y20A) *
                    std::exp(-rA2 / sig2);

                amrex::Real delayed_kick = kickA;
                if (metric_type != 1)
                {
                    const amrex::Real dxB = xloc - centerB[0];
                    const amrex::Real dyB = yloc - centerB[1];
                    const amrex::Real dzB = zloc - centerB[2];
                    const amrex::Real rB2 = dxB * dxB + dyB * dyB + dzB * dzB;
                    const amrex::Real rB2_safe =
                        amrex::max(rB2, amrex::Real(1.0e-24));
                    const amrex::Real y20B =
                        3.0 * dzB * dzB / rB2_safe - 1.0;
                    delayed_kick +=
                        (delayed_monopole + delayed_quadrupole * y20B) *
                        std::exp(-rB2 / sig2);
                }

                rhs_arrs[box_no](i, j, k, c_K) +=
                    temporal_weight * delayed_kick;
            });
    }

    amrex::Gpu::streamSynchronize();
}

void WormholeLevel::specificUpdateODE(amrex::MultiFab &a_soln)
{
    const auto &soln_arrs = a_soln.arrays();
    TraceARemoval trace_A_removal;
    amrex::ParallelFor(a_soln, amrex::IntVect(0),
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       {
                           trace_A_removal(i, j, k, soln_arrs[box_no]);
                       });

    amrex::Gpu::streamSynchronize();
}

void WormholeLevel::pre_tag_cells()
{
    amrex::MultiFab &state_new = get_new_data(state_index);
    const auto cur_time        = get_state_data(state_index).curTime();
    FillPatch(*this, state_new, 2, cur_time, state_index, c_chi, 1);
}

void WormholeLevel::tag_cells(amrex::TagBoxArray &a_tag_box_array,
                              amrex::Real a_regrid_threshold)
{
    BL_PROFILE("WormholeLevel::tag_cells()");
    amrex::MultiFab &state_new = get_new_data(state_index);
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
        const amrex::Real time         = get_state_data(state_index).curTime();
        const amrex::Real dt           = parent->dtLevel(0);
        const amrex::Real restart_time = get_gramr_ptr()->get_restart_time();
        const bool first_step          = (time == 0.0);

        // Fill ghosts for constraint calculation
        amrex::MultiFab &state_new = get_new_data(state_index);
        FillPatch(*this, state_new, 2, time, state_index, 0, state_new.nComp());

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

    // --- Collapse / strong-field diagnostics (small-data output) ---
    // These are cheap global reductions that can be recorded even when plotfiles
    // are deleted on-the-fly.
    if (Level() == 0)
    {
        const amrex::Real time         = get_state_data(state_index).curTime();
        const amrex::Real dt           = parent->dtLevel(0);
        const amrex::Real restart_time = get_gramr_ptr()->get_restart_time();
        const bool first_step          = (time == 0.0);

        amrex::MultiFab &state_new = get_new_data(state_index);
        FillPatch(*this, state_new, 2, time, state_index, 0, state_new.nComp());

        amrex::ReduceOps<amrex::ReduceOpMin, amrex::ReduceOpMin,
                         amrex::ReduceOpMax, amrex::ReduceOpMax>
            reduce_ops;
        amrex::ReduceData<amrex::Real, amrex::Real, amrex::Real, amrex::Real>
            reduce_data(reduce_ops);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        const auto dx_arr  = Geom().CellSizeArray();
        const auto prob_lo = Geom().ProbLoArray();

        for (amrex::MFIter mfi(state_new, amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            const auto arr       = state_new.const_array(mfi);
            reduce_ops.eval(
                bx, reduce_data,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple
                {
                    const amrex::Real lapse = arr(i, j, k, c_lapse);
                    const amrex::Real chi   = arr(i, j, k, c_chi);
                    const amrex::Real K     = arr(i, j, k, c_K);

                    const amrex::Real x =
                        prob_lo[0] + (amrex::Real(i) + 0.5) * dx_arr[0];
                    const amrex::Real y =
                        prob_lo[1] + (amrex::Real(j) + 0.5) * dx_arr[1];
                    const amrex::Real z =
                        prob_lo[2] + (amrex::Real(k) + 0.5) * dx_arr[2];
                    const amrex::Real r2 = x * x + y * y + z * z;
                    const amrex::Real r  = std::sqrt(r2);

                    amrex::Real ah_radius = 0.0;
                    if (r > 1e-6) // avoid division by zero at the throat center
                    {
                        const amrex::Real A11 = arr(i, j, k, c_A11);
                        const amrex::Real A22 = arr(i, j, k, c_A22);
                        const amrex::Real A33 = arr(i, j, k, c_A33);
                        const amrex::Real A12 = arr(i, j, k, c_A12);
                        const amrex::Real A13 = arr(i, j, k, c_A13);
                        const amrex::Real A23 = arr(i, j, k, c_A23);

                        const amrex::Real Arr =
                            (A11 * x * x + A22 * y * y + A33 * z * z +
                             2.0 * A12 * x * y + 2.0 * A13 * x * z +
                             2.0 * A23 * y * z) /
                            r2;

                        // Spherical trapped-surface proxy:
                        // theta_+ = 2 sqrt(chi) / r + 2/3 K - chi * A_rr
                        const amrex::Real theta_plus =
                            2.0 * std::sqrt(chi) / r + (2.0 / 3.0) * K -
                            chi * Arr;

                        if (theta_plus <= 0.0)
                        {
                            ah_radius = r;
                        }
                    }

                    return {lapse, chi, amrex::Math::abs(K), ah_radius};
                });
        }

        const auto reduce_vals = reduce_data.value();
        amrex::Real min_lapse  = amrex::get<0>(reduce_vals);
        amrex::Real min_chi    = amrex::get<1>(reduce_vals);
        amrex::Real max_abs_K  = amrex::get<2>(reduce_vals);
        amrex::Real max_ah_r   = amrex::get<3>(reduce_vals);
        amrex::ParallelDescriptor::ReduceRealMin(min_lapse);
        amrex::ParallelDescriptor::ReduceRealMin(min_chi);
        amrex::ParallelDescriptor::ReduceRealMax(max_abs_K);
        amrex::ParallelDescriptor::ReduceRealMax(max_ah_r);

        // Location of the global minimum lapse (average over ties).
        // This helps confirm the collapse localizes at the throat/center rather
        // than being a boundary artifact.
        const amrex::Real tol =
            amrex::max(amrex::Real(1.0e-14), amrex::Real(1.0e-12) * amrex::Math::abs(min_lapse));

        amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpSum,
                         amrex::ReduceOpSum, amrex::ReduceOpSum>
            reduce_ops_loc;
        amrex::ReduceData<amrex::Real, amrex::Real, amrex::Real, amrex::Real>
            reduce_data_loc(reduce_ops_loc);
        using ReduceTupleLoc = typename decltype(reduce_data_loc)::Type;

        for (amrex::MFIter mfi(state_new, amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            const auto arr       = state_new.const_array(mfi);
            reduce_ops_loc.eval(
                bx, reduce_data_loc,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTupleLoc
                {
                    const amrex::Real lapse = arr(i, j, k, c_lapse);
                    const bool is_min       = (amrex::Math::abs(lapse - min_lapse) <= tol);
                    if (!is_min)
                    {
                        return {0.0, 0.0, 0.0, 0.0};
                    }
                    const amrex::Real x =
                        prob_lo[0] + (amrex::Real(i) + 0.5) * dx[0];
                    const amrex::Real y =
                        prob_lo[1] + (amrex::Real(j) + 0.5) * dx[1];
                    const amrex::Real z =
                        prob_lo[2] + (amrex::Real(k) + 0.5) * dx[2];
                    return {x, y, z, 1.0};
                });
        }

        auto [sum_x, sum_y, sum_z, count] = reduce_data_loc.value();
        amrex::ParallelDescriptor::ReduceRealSum(sum_x);
        amrex::ParallelDescriptor::ReduceRealSum(sum_y);
        amrex::ParallelDescriptor::ReduceRealSum(sum_z);
        amrex::ParallelDescriptor::ReduceRealSum(count);

        const amrex::Real min_lapse_x = (count > 0.0) ? (sum_x / count) : 0.0;
        const amrex::Real min_lapse_y = (count > 0.0) ? (sum_y / count) : 0.0;
        const amrex::Real min_lapse_z = (count > 0.0) ? (sum_z / count) : 0.0;

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

        const std::string prefix = out_dir + "collapse_diagnostics";
        SmallDataIO diag_file(prefix, dt, time, restart_time, SmallDataIO::APPEND,
                              first_step);
        diag_file.remove_duplicate_time_data();
        if (first_step)
        {
            diag_file.write_header_line(
                {"min_lapse", "min_chi", "max_abs_K", "min_lapse_x",
                 "min_lapse_y", "min_lapse_z", "max_ah_r"});
        }
        diag_file.write_time_data_line({static_cast<double>(min_lapse),
                                        static_cast<double>(min_chi),
                                        static_cast<double>(max_abs_K),
                                        static_cast<double>(min_lapse_x),
                                        static_cast<double>(min_lapse_y),
                                        static_cast<double>(min_lapse_z),
                                        static_cast<double>(max_ah_r)});
    }
}