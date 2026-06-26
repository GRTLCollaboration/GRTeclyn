#ifndef RADIALRECIPE_CONSTRAINT_NORMS_HPP_
#define RADIALRECIPE_CONSTRAINT_NORMS_HPP_

#include "ComplexScalarField.hpp"
#include "ConstraintsWithMatter.hpp"
#include "ExoticScalarField.hpp"
#include "GRTresnaIndependentScalars.hpp"
#include "RadialRecipeLevel.hpp"
#include "RadialRecipeMatterConstraints.hpp"
#include "RadialRecipeMatterDispatch.hpp"
#include "ScalarField.hpp"
#include "SetupFunctions.hpp"

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

    if (RadialRecipeMatter::uses_independent_scalars(sim_params))
    {
        GRTresnaIndependentScalars matter(
            sim_params.recipe_num_scalar_fields,
            sim_params.recipe_scalar_field_signs, sim_params.recipe_scalar_mass,
            sim_params.recipe_scalar_lambda);
        fill_matter_constraints(cst, state_new, matter, dx[0],
                                sim_params.recipe_params.grid_center, time);
    }
    else if (RadialRecipeMatter::uses_complex_scalar(sim_params))
    {
        ComplexScalarField matter(sim_params.recipe_scalar_mass,
                                  sim_params.recipe_scalar_lambda,
                                  sim_params.recipe_scalar_sign);
        fill_matter_constraints(cst, state_new, matter, dx[0],
                                sim_params.recipe_params.grid_center, time);
    }
    else if (RadialRecipeMatter::uses_bicomplex_scalar(sim_params))
    {
        BiComplexScalarField matter(sim_params.recipe_scalar_mass,
                                    sim_params.recipe_scalar_lambda);
        fill_matter_constraints(cst, state_new, matter, dx[0],
                                sim_params.recipe_params.grid_center, time);
    }
    else if (sim_params.recipe_exotic_matter)
    {
        ExoticScalarField<DefaultPotential> exotic_scalar(
            DefaultPotential(), sim_params.recipe_support_strength);
        fill_matter_constraints(cst, state_new, exotic_scalar, dx[0],
                                sim_params.recipe_params.grid_center, time);
    }
    else
    {
        ScalarField<DefaultPotential> scalar_field;
        fill_matter_constraints(cst, state_new, scalar_field, dx[0],
                                sim_params.recipe_params.grid_center, time);
    }

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
