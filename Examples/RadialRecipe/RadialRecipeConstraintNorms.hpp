#ifndef RADIALRECIPE_CONSTRAINT_NORMS_HPP_
#define RADIALRECIPE_CONSTRAINT_NORMS_HPP_

#include "RadialRecipeLevel.hpp"
#include "RadialRecipeMatterDispatch.hpp"
#include "SetupFunctions.hpp"
#include "StateTypes.hpp"

#include <AMReX_MultiFab.H>
#include <AMReX_Reduce.H>

struct RadialRecipeConstraintNorms
{
    double L2_Ham{0.0};
    double L2_Mom{0.0};
};

inline RadialRecipeConstraintNorms
compute_radial_recipe_constraint_norms(RadialRecipeLevel &level)
{
    RadialRecipeConstraintNorms out;
    const auto &sim_params = RadialRecipeLevel::simParams();

    const amrex::Real time = level.get_state_data(state_index).curTime();
    amrex::MultiFab &state_new = level.get_new_data(state_index);
    amrex::AmrLevel::FillPatch(level, state_new, 2, time, state_index, 0,
                               state_new.nComp());

    amrex::MultiFab cst(state_new.boxArray(), state_new.DistributionMap(), 4, 0);
    cst.setVal(0.0);
    const auto dx = level.Geom().CellSizeArray();

    RadialRecipeMatter::fill_active_constraints(
        cst, state_new, sim_params, dx[0], sim_params.recipe_params.grid_center,
        time, /*with_abs_terms=*/false);

    const amrex::Real cell_vol = dx[0] * dx[1] * dx[2];

    amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum>
        reduce_ops;
    amrex::ReduceData<amrex::Real, amrex::Real, amrex::Real> reduce_data(
        reduce_ops);

    for (amrex::MFIter mfi(cst, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.validbox();
        const auto arr         = cst.const_array(mfi);
        reduce_ops.eval(
            bx, reduce_data,
            [=] AMREX_GPU_DEVICE(int i, int j, int k)
                -> amrex::GpuTuple<amrex::Real, amrex::Real, amrex::Real>
            {
                const amrex::Real ham  = arr(i, j, k, 0);
                const amrex::Real m1   = arr(i, j, k, 1);
                const amrex::Real m2   = arr(i, j, k, 2);
                const amrex::Real m3   = arr(i, j, k, 3);
                const amrex::Real mom2 = (m1 * m1 + m2 * m2 + m3 * m3);
                return {ham * ham * cell_vol, mom2 * cell_vol, cell_vol};
            });
    }

    auto [sum_ham2, sum_mom2, sum_vol] = reduce_data.value();
    amrex::ParallelDescriptor::ReduceRealSum(sum_ham2);
    amrex::ParallelDescriptor::ReduceRealSum(sum_mom2);
    amrex::ParallelDescriptor::ReduceRealSum(sum_vol);

    out.L2_Ham = (sum_vol > 0.0) ? std::sqrt(sum_ham2 / sum_vol) : 0.0;
    out.L2_Mom = (sum_vol > 0.0) ? std::sqrt(sum_mom2 / sum_vol) : 0.0;
    return out;
}

#endif /* RADIALRECIPE_CONSTRAINT_NORMS_HPP_ */
