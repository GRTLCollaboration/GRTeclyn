#ifndef RL_OBSERVATION_COLLECTOR_HPP_
#define RL_OBSERVATION_COLLECTOR_HPP_

#include "GRAMR.hpp"
#include "RadialRecipeConstraintNorms.hpp"
#include "RadialRecipeLevel.hpp"
#include "SetupFunctions.hpp"
#include "StateVariables.hpp"

#include <AMReX_Reduce.H>
#include <array>

inline std::array<double, 6> collect_rl_observations(GRAMR &amr,
                                                     double l2_ham, double l2_mom)
{
    std::array<double, 6> obs{};

    const int finest_lev = amr.finestLevel();
    auto &fine_level     = dynamic_cast<RadialRecipeLevel &>(amr.getLevel(finest_lev));
    amrex::MultiFab &state_fine =
        fine_level.get_new_data(RadialRecipeLevel::state_index);
    const amrex::Real time =
        fine_level.get_state_data(RadialRecipeLevel::state_index).curTime();

    FillPatch(fine_level, state_fine, 2, time, RadialRecipeLevel::state_index, 0,
              state_fine.nComp());

    amrex::ReduceOps<amrex::ReduceOpMin, amrex::ReduceOpMin, amrex::ReduceOpMax>
        reduce_ops;
    amrex::ReduceData<amrex::Real, amrex::Real, amrex::Real> reduce_data(
        reduce_ops);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    for (amrex::MFIter mfi(state_fine, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi)
    {
        const amrex::Box &bx = mfi.validbox();
        const auto arr       = state_fine.const_array(mfi);
        reduce_ops.eval(
            bx, reduce_data,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple
            {
                const amrex::Real chi   = arr(i, j, k, c_chi);
                const amrex::Real lapse = arr(i, j, k, c_lapse);
                const amrex::Real K     = arr(i, j, k, c_K);
                return {chi, lapse, std::abs(K)};
            });
    }

    auto [min_chi, min_lapse, max_abs_K] = reduce_data.value();
    amrex::ParallelDescriptor::ReduceRealMin(min_chi);
    amrex::ParallelDescriptor::ReduceRealMin(min_lapse);
    amrex::ParallelDescriptor::ReduceRealMax(max_abs_K);

    obs[0] = static_cast<double>(min_chi);
    obs[1] = static_cast<double>(min_lapse);
    obs[2] = static_cast<double>(max_abs_K);
    obs[3] = l2_ham;
    obs[4] = l2_mom;
    obs[5] = static_cast<double>(time);
    return obs;
}

inline double compute_l2_hamiltonian_norm(GRAMR &amr)
{
    auto &level0 = dynamic_cast<RadialRecipeLevel &>(amr.getLevel(0));
    return compute_radial_recipe_constraint_norms(level0).L2_Ham;
}

#endif /* RL_OBSERVATION_COLLECTOR_HPP_ */
