#include "RadialRecipeInitialData.hpp"
#include "RadialRecipeLevel.hpp"
#include "CCZ4RHSWithMatter.hpp"
#include "ChiTagger.hpp"
#include "ConstraintsWithMatter.hpp"
#include "GRParmParse.hpp"
#include "PositiveChiAndLapse.hpp"
#include "ScalarField.hpp"
#include "SmallDataIO.hpp"
#include "TraceARemoval.hpp"
#include "Weyl4WithMatter.hpp"

#include <AMReX_Reduce.H>
#include <AMReX_Utility.H>
#include <cmath>

void RadialRecipeLevel::variableSetUp()
{
    BL_PROFILE("RadialRecipeLevel::variableSetUp()");
    stateVariableSetUp();

    ConstraintsWithMatter<ScalarField<DefaultPotential>>::set_up(state_index);
    Weyl4WithMatter<ScalarField<DefaultPotential>>::set_up(state_index);
}

void RadialRecipeLevel::specificAdvance()
{
    amrex::MultiFab &S_new = get_new_data(state_index);
    const auto &arrs       = S_new.arrays();
    TraceARemoval trace_A_removal;
    PositiveChiAndLapse positive_chi_lapse(simParams().min_chi,
                                           simParams().min_lapse);

    amrex::ParallelFor(S_new,
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       {
                           trace_A_removal(i, j, k, arrs[box_no]);
                           positive_chi_lapse(i, j, k, arrs[box_no]);
                       });
}

void RadialRecipeLevel::initData()
{
    BL_PROFILE("RadialRecipeLevel::initData");

    RadialRecipeInitialData recipe(simParams().recipe_params, Geom().CellSize(0));

    amrex::MultiFab &state = get_new_data(state_index);
    const auto &arrs       = state.arrays();

    amrex::ParallelFor(state, state.nGrowVect(),
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       {
                           amrex::CellData<amrex::Real> cell =
                               arrs[box_no].cellData(i, j, k);
                           for (int n = 0; n < cell.nComp(); ++n)
                           {
                               cell[n] = 0.;
                           }
                           recipe.compute(i, j, k, arrs[box_no]);
                       });

    amrex::Gpu::streamSynchronize();
}

void RadialRecipeLevel::specificEvalRHS(amrex::MultiFab &a_soln,
                                        amrex::MultiFab &a_rhs,
                                        const double a_time)
{
    BL_PROFILE("RadialRecipeLevel::specificEvalRHS()");
    const int soln_ghosts = a_soln.nGrowVect()[0];
    if (soln_ghosts > 0)
    {
        FillPatch(*this, a_soln, soln_ghosts, a_time, state_index, 0,
                  a_soln.nComp());
    }
    const auto &soln_arrs   = a_soln.arrays();
    const auto &soln_c_arrs = a_soln.const_arrays();
    const auto &rhs_arrs    = a_rhs.arrays();
    TraceARemoval trace_A_removal;
    PositiveChiAndLapse positive_chi_lapse(simParams().min_chi,
                                           simParams().min_lapse);

    amrex::ParallelFor(a_soln, a_soln.nGrowVect(),
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       {
                           trace_A_removal(i, j, k, soln_arrs[box_no]);
                           positive_chi_lapse(i, j, k, soln_arrs[box_no]);
                       });

    ScalarField<DefaultPotential> scalar_field;
    CCZ4RHSWithMatter<ScalarField<DefaultPotential>, MovingPunctureGaugeWithMatter,
                      FourthOrderDerivatives>
        ccz4rhs(scalar_field, simParams().ccz4_params, Geom().CellSize(0),
                simParams().sigma, simParams().formulation, 1.0,
                simParams().recipe_params.grid_center, a_time);

    amrex::ParallelFor(
        a_rhs, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
        { ccz4rhs(i, j, k, rhs_arrs[box_no], soln_c_arrs[box_no]); });

    amrex::Gpu::streamSynchronize();
}

void RadialRecipeLevel::specificUpdateODE(amrex::MultiFab &a_soln)
{
    const auto &soln_arrs = a_soln.arrays();
    TraceARemoval trace_A_removal;
    PositiveChiAndLapse positive_chi_lapse(simParams().min_chi,
                                           simParams().min_lapse);
    amrex::ParallelFor(a_soln, amrex::IntVect(0),
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       {
                           trace_A_removal(i, j, k, soln_arrs[box_no]);
                           positive_chi_lapse(i, j, k, soln_arrs[box_no]);
                       });

    amrex::Gpu::streamSynchronize();
}

void RadialRecipeLevel::pre_tag_cells()
{
    amrex::MultiFab &state_new = get_new_data(state_index);
    const auto cur_time        = get_state_data(state_index).curTime();
    FillPatch(*this, state_new, 2, cur_time, state_index, c_chi, 1);
}

void RadialRecipeLevel::tag_cells(amrex::TagBoxArray &a_tag_box_array,
                                  amrex::Real a_regrid_threshold)
{
    BL_PROFILE("RadialRecipeLevel::tag_cells()");
    amrex::MultiFab &state_new = get_new_data(state_index);
    const auto &tag_arrs       = a_tag_box_array.arrays();
    const auto &state_new_arrs = state_new.const_arrays();

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

void RadialRecipeLevel::specificPostTimeStep()
{
    BL_PROFILE("RadialRecipeLevel::specificPostTimeStep");

    if (simParams().calculate_constraint_norms && Level() == 0)
    {
        const amrex::Real time         = get_state_data(state_index).curTime();
        const amrex::Real dt           = parent->dtLevel(0);
        const amrex::Real restart_time = get_gramr_ptr()->get_restart_time();
        const bool first_step          = (time == 0.0);

        amrex::MultiFab &state_new = get_new_data(state_index);
        FillPatch(*this, state_new, 2, time, state_index, 0, state_new.nComp());

        amrex::MultiFab cst(state_new.boxArray(), state_new.DistributionMap(), 4,
                            0);
        cst.setVal(0.0);
        ScalarField<DefaultPotential> scalar_field;
        const auto dx = Geom().CellSizeArray();
        ConstraintsWithMatter<ScalarField<DefaultPotential>> my_constraints(
            scalar_field, dx[0], 1.0, 0, Interval(1, 3),
            simParams().recipe_params.grid_center, time);

        for (amrex::MFIter mfi(cst, amrex::TilingIfNotGPU()); mfi.isValid();
             ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            const auto arr       = cst.array(mfi);
            const auto src_arr   = state_new.const_array(mfi);

            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int ix, int iy, int iz) noexcept
                { my_constraints(ix, iy, iz, arr, src_arr); });
        }

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

        amrex::MultiFab cst_vac(state_new.boxArray(),
                                state_new.DistributionMap(), 4, 0);
        cst_vac.setVal(0.0);
        Constraints vacuum_constraints(dx[0], 0, Interval(1, 3));
        for (amrex::MFIter mfi(cst_vac, amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            const auto arr       = cst_vac.array(mfi);
            const auto src_arr   = state_new.const_array(mfi);
            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int ix, int iy, int iz) noexcept
                { vacuum_constraints(ix, iy, iz, arr, src_arr); });
        }

        constexpr amrex::Real inv_16pi = 1.0 / (16.0 * M_PI);

        amrex::ReduceOps<amrex::ReduceOpMin, amrex::ReduceOpMax,
                         amrex::ReduceOpSum, amrex::ReduceOpSum>
            rho_reduce_ops;
        amrex::ReduceData<amrex::Real, amrex::Real, amrex::Real, amrex::Real>
            rho_reduce_data(rho_reduce_ops);
        using RhoReduceTuple = typename decltype(rho_reduce_data)::Type;

        for (amrex::MFIter mfi(cst_vac, amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            const auto arr       = cst_vac.const_array(mfi);
            rho_reduce_ops.eval(
                bx, rho_reduce_data,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) -> RhoReduceTuple
                {
                    const amrex::Real ham_vac = arr(i, j, k, 0);
                    const amrex::Real rho_req = ham_vac * inv_16pi;
                    const amrex::Real neg_rho =
                        (rho_req < 0.0) ? (-rho_req * cell_vol) : 0.0;
                    return {rho_req, rho_req, neg_rho, cell_vol};
                });
        }

        auto rho_vals = rho_reduce_data.value();
        amrex::Real min_rho_req     = amrex::get<0>(rho_vals);
        amrex::Real max_rho_req     = amrex::get<1>(rho_vals);
        amrex::Real sum_neg_rho     = amrex::get<2>(rho_vals);
        amrex::Real sum_rho_vol     = amrex::get<3>(rho_vals);
        amrex::ParallelDescriptor::ReduceRealMin(min_rho_req);
        amrex::ParallelDescriptor::ReduceRealMax(max_rho_req);
        amrex::ParallelDescriptor::ReduceRealSum(sum_neg_rho);
        amrex::ParallelDescriptor::ReduceRealSum(sum_rho_vol);

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
            constraints_file.write_header_line(
                {"L2_Ham", "L2_Mom", "min_rho_req", "max_rho_req",
                 "integral_neg_rho"});
        }
        constraints_file.write_time_data_line(
            {L2_Ham, L2_Mom, static_cast<double>(min_rho_req),
             static_cast<double>(max_rho_req),
             static_cast<double>(sum_neg_rho)});
    }

    if (Level() == 0)
    {
        const amrex::Real time         = get_state_data(state_index).curTime();
        const amrex::Real dt           = parent->dtLevel(0);
        const amrex::Real restart_time = get_gramr_ptr()->get_restart_time();
        const bool first_step          = (time == 0.0);

        const int finest_lev = parent->finestLevel();
        auto &fine_level     = parent->getLevel(finest_lev);
        amrex::MultiFab &state_fine = fine_level.get_new_data(state_index);
        const auto &fine_geom       = parent->Geom(finest_lev);

        FillPatch(fine_level, state_fine, 2, time, state_index, 0,
                  state_fine.nComp());

        {
            const auto &arrs = state_fine.arrays();
            TraceARemoval trace_A_removal;
            PositiveChiAndLapse positive_chi_lapse(simParams().min_chi,
                                                   simParams().min_lapse);
            amrex::ParallelFor(state_fine, amrex::IntVect(0),
                               [=] AMREX_GPU_DEVICE(int box_no, int i, int j,
                                                    int k)
                               {
                                   trace_A_removal(i, j, k, arrs[box_no]);
                                   positive_chi_lapse(i, j, k, arrs[box_no]);
                               });
            amrex::Gpu::streamSynchronize();
        }

        amrex::ReduceOps<amrex::ReduceOpMin, amrex::ReduceOpMin,
                         amrex::ReduceOpMax, amrex::ReduceOpMax,
                         amrex::ReduceOpMin,
                         amrex::ReduceOpMin, amrex::ReduceOpMax,
                         amrex::ReduceOpMin, amrex::ReduceOpMax>
            reduce_ops;
        amrex::ReduceData<amrex::Real, amrex::Real, amrex::Real, amrex::Real,
                          amrex::Real,
                          amrex::Real, amrex::Real,
                          amrex::Real, amrex::Real>
            reduce_data(reduce_ops);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        const auto prob_lo = fine_geom.ProbLoArray();
        const auto dx_arr = fine_geom.CellSizeArray();

        for (amrex::MFIter mfi(state_fine, amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            const auto arr       = state_fine.const_array(mfi);
            reduce_ops.eval(
                bx, reduce_data,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple
                {
                    const amrex::Real lapse = arr(i, j, k, c_lapse);
                    const amrex::Real chi   = arr(i, j, k, c_chi);
                    const amrex::Real K     = arr(i, j, k, c_K);

                    const amrex::Real x = prob_lo[0] + (amrex::Real(i) + 0.5) * dx_arr[0];
                    const amrex::Real y = prob_lo[1] + (amrex::Real(j) + 0.5) * dx_arr[1];
                    const amrex::Real z = prob_lo[2] + (amrex::Real(k) + 0.5) * dx_arr[2];
                    const amrex::Real r2 = x*x + y*y + z*z;
                    const amrex::Real r = std::sqrt(r2);

                    amrex::Real ah_radius = 0.0;
                    amrex::Real theta_plus_min_proxy = 1.0e30;
                    if (r > 1e-6)
                    {
                        const amrex::Real A11 = arr(i, j, k, c_A11);
                        const amrex::Real A22 = arr(i, j, k, c_A22);
                        const amrex::Real A33 = arr(i, j, k, c_A33);
                        const amrex::Real A12 = arr(i, j, k, c_A12);
                        const amrex::Real A13 = arr(i, j, k, c_A13);
                        const amrex::Real A23 = arr(i, j, k, c_A23);

                        const amrex::Real Arr = (A11*x*x + A22*y*y + A33*z*z + 2.0*A12*x*y + 2.0*A13*x*z + 2.0*A23*y*z) / r2;

                        const amrex::Real dx_chi =
                            (arr(i + 1, j, k, c_chi) - arr(i - 1, j, k, c_chi)) /
                            (2.0 * dx_arr[0]);
                        const amrex::Real dy_chi =
                            (arr(i, j + 1, k, c_chi) - arr(i, j - 1, k, c_chi)) /
                            (2.0 * dx_arr[1]);
                        const amrex::Real dz_chi =
                            (arr(i, j, k + 1, c_chi) - arr(i, j, k - 1, c_chi)) /
                            (2.0 * dx_arr[2]);
                        const amrex::Real dchi_dr =
                            (x * dx_chi + y * dy_chi + z * dz_chi) / r;
                        const amrex::Real sqrt_chi =
                            std::sqrt(amrex::max(chi, amrex::Real(1.0e-20)));

                        const amrex::Real theta_plus =
                            2.0 * sqrt_chi / r - dchi_dr / sqrt_chi + Arr -
                            (2.0 / 3.0) * K;
                        theta_plus_min_proxy = theta_plus;

                        if (theta_plus <= 0.0)
                        {
                            ah_radius = r;
                        }
                    }

                    const amrex::Real sf_phi = arr(i, j, k, c_phi);
                    const amrex::Real sf_Pi  = arr(i, j, k, c_Pi);

                    return {lapse, chi, amrex::Math::abs(K), ah_radius, theta_plus_min_proxy,
                            sf_phi, sf_phi, sf_Pi, sf_Pi};
                });
        }

        const auto reduce_vals = reduce_data.value();
        amrex::Real min_lapse  = amrex::get<0>(reduce_vals);
        amrex::Real min_chi    = amrex::get<1>(reduce_vals);
        amrex::Real max_abs_K  = amrex::get<2>(reduce_vals);
        amrex::Real max_ah_r   = amrex::get<3>(reduce_vals);
        amrex::Real min_theta_plus = amrex::get<4>(reduce_vals);
        amrex::Real min_phi    = amrex::get<5>(reduce_vals);
        amrex::Real max_phi    = amrex::get<6>(reduce_vals);
        amrex::Real min_Pi     = amrex::get<7>(reduce_vals);
        amrex::Real max_Pi     = amrex::get<8>(reduce_vals);
        amrex::ParallelDescriptor::ReduceRealMin(min_lapse);
        amrex::ParallelDescriptor::ReduceRealMin(min_chi);
        amrex::ParallelDescriptor::ReduceRealMax(max_abs_K);
        amrex::ParallelDescriptor::ReduceRealMax(max_ah_r);
        amrex::ParallelDescriptor::ReduceRealMin(min_theta_plus);
        amrex::ParallelDescriptor::ReduceRealMin(min_phi);
        amrex::ParallelDescriptor::ReduceRealMax(max_phi);
        amrex::ParallelDescriptor::ReduceRealMin(min_Pi);
        amrex::ParallelDescriptor::ReduceRealMax(max_Pi);

        const amrex::Real tol =
            amrex::max(amrex::Real(1.0e-14), amrex::Real(1.0e-12) * amrex::Math::abs(min_lapse));

        amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpSum,
                         amrex::ReduceOpSum, amrex::ReduceOpSum>
            reduce_ops_loc;
        amrex::ReduceData<amrex::Real, amrex::Real, amrex::Real, amrex::Real>
            reduce_data_loc(reduce_ops_loc);
        using ReduceTupleLoc = typename decltype(reduce_data_loc)::Type;

        for (amrex::MFIter mfi(state_fine, amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            const auto arr       = state_fine.const_array(mfi);
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
                        prob_lo[0] + (amrex::Real(i) + 0.5) * dx_arr[0];
                    const amrex::Real y =
                        prob_lo[1] + (amrex::Real(j) + 0.5) * dx_arr[1];
                    const amrex::Real z =
                        prob_lo[2] + (amrex::Real(k) + 0.5) * dx_arr[2];
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

        const amrex::Real tol_theta = amrex::max(
            amrex::Real(1.0e-12),
            amrex::Real(1.0e-8) * amrex::Math::abs(min_theta_plus));

        amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpSum>
            reduce_ops_theta_loc;
        amrex::ReduceData<amrex::Real, amrex::Real> reduce_data_theta_loc(
            reduce_ops_theta_loc);
        using ReduceTupleThetaLoc = typename decltype(reduce_data_theta_loc)::Type;

        for (amrex::MFIter mfi(state_fine, amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            const auto arr       = state_fine.const_array(mfi);
            reduce_ops_theta_loc.eval(
                bx, reduce_data_theta_loc,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTupleThetaLoc
                {
                    const amrex::Real chi = arr(i, j, k, c_chi);
                    const amrex::Real K   = arr(i, j, k, c_K);

                    const amrex::Real x =
                        prob_lo[0] + (amrex::Real(i) + 0.5) * dx_arr[0];
                    const amrex::Real y =
                        prob_lo[1] + (amrex::Real(j) + 0.5) * dx_arr[1];
                    const amrex::Real z =
                        prob_lo[2] + (amrex::Real(k) + 0.5) * dx_arr[2];
                    const amrex::Real r2 = x * x + y * y + z * z;
                    const amrex::Real r  = std::sqrt(r2);

                    if (r <= 1e-6)
                    {
                        return {0.0, 0.0};
                    }

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

                    const amrex::Real dx_chi =
                        (arr(i + 1, j, k, c_chi) - arr(i - 1, j, k, c_chi)) /
                        (2.0 * dx_arr[0]);
                    const amrex::Real dy_chi =
                        (arr(i, j + 1, k, c_chi) - arr(i, j - 1, k, c_chi)) /
                        (2.0 * dx_arr[1]);
                    const amrex::Real dz_chi =
                        (arr(i, j, k + 1, c_chi) - arr(i, j, k - 1, c_chi)) /
                        (2.0 * dx_arr[2]);
                    const amrex::Real dchi_dr =
                        (x * dx_chi + y * dy_chi + z * dz_chi) / r;
                    const amrex::Real sqrt_chi =
                        std::sqrt(amrex::max(chi, amrex::Real(1.0e-20)));

                    const amrex::Real theta_plus =
                        2.0 * sqrt_chi / r - dchi_dr / sqrt_chi + Arr -
                        (2.0 / 3.0) * K;

                    const bool is_min =
                        (amrex::Math::abs(theta_plus - min_theta_plus) <=
                         tol_theta);
                    if (!is_min)
                    {
                        return {0.0, 0.0};
                    }
                    return {r, 1.0};
                });
        }

        auto [sum_r_theta, count_r_theta] = reduce_data_theta_loc.value();
        amrex::ParallelDescriptor::ReduceRealSum(sum_r_theta);
        amrex::ParallelDescriptor::ReduceRealSum(count_r_theta);
        const amrex::Real r_at_min_theta_plus =
            (count_r_theta > 0.0) ? (sum_r_theta / count_r_theta) : 0.0;

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
                 "min_lapse_y", "min_lapse_z", "max_ah_r", "min_theta_plus",
                 "r_at_min_theta_plus",
                 "min_phi", "max_phi", "min_Pi", "max_Pi"});
        }
        diag_file.write_time_data_line({static_cast<double>(min_lapse),
                                        static_cast<double>(min_chi),
                                        static_cast<double>(max_abs_K),
                                        static_cast<double>(min_lapse_x),
                                        static_cast<double>(min_lapse_y),
                                        static_cast<double>(min_lapse_z),
                                        static_cast<double>(max_ah_r),
                                        static_cast<double>(min_theta_plus),
                                        static_cast<double>(r_at_min_theta_plus),
                                        static_cast<double>(min_phi),
                                        static_cast<double>(max_phi),
                                        static_cast<double>(min_Pi),
                                        static_cast<double>(max_Pi)});
    }
}
