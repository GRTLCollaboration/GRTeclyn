#ifndef RL_OBSERVATION_COLLECTOR_HPP_
#define RL_OBSERVATION_COLLECTOR_HPP_

#include "GRAMR.hpp"
#include "RLLumpState.hpp"
#include "RLLumpTracker.hpp"
#include "RLMatterPumpParams.hpp" // RLRuntime::g_lump_state
#include "RadialRecipeConstraintNorms.hpp"
#include "RadialRecipeLevel.hpp"
#include "SetupFunctions.hpp"
#include "StateVariables.hpp"

#include <AMReX_Reduce.H>
#include <array>
#include <vector>

//! Build the RL observation vector: 6 global survival/constraint scalars
//! followed by ``RL_LUMP_OBS_STRIDE`` per-lump entries for each tracked lump
//! (so the agent sees where every lump is, its size, mass and local collapse
//! state in 3-D).  Runs the lump tracker as a side effect, updating the
//! persistent ``RLRuntime::g_lump_state`` (centres also drive the pump).
//!
//! Layout: [min_chi, min_lapse, max_abs_K, L2_Ham, L2_Mom, time,
//!          {x,y,z,size,mass,peak,min_lapse,min_chi} * num_lumps]
inline std::vector<double> collect_rl_observations(
    GRAMR &amr, double l2_ham, double l2_mom, int num_lumps, bool complex_field,
    double ball_radius, const std::array<double, AMREX_SPACEDIM> &center)
{
    const int finest_lev = amr.finestLevel();
    auto &fine_level = dynamic_cast<RadialRecipeLevel &>(amr.getLevel(finest_lev));
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

    // Per-lump 3-D state (also moves the pump spotlights to follow the lumps).
    track_rl_lumps(amr, num_lumps, complex_field, ball_radius, center,
                   RLRuntime::g_lump_state);

    if (num_lumps < 0)
        num_lumps = 0;
    if (num_lumps > RL_MAX_LUMPS)
        num_lumps = RL_MAX_LUMPS;

    std::vector<double> obs;
    obs.reserve(static_cast<std::size_t>(6 + RL_LUMP_OBS_STRIDE * num_lumps));
    obs.push_back(static_cast<double>(min_chi));
    obs.push_back(static_cast<double>(min_lapse));
    obs.push_back(static_cast<double>(max_abs_K));
    obs.push_back(l2_ham);
    obs.push_back(l2_mom);
    obs.push_back(static_cast<double>(time));
    for (int s = 0; s < num_lumps; ++s)
    {
        const RLLumpState &L = RLRuntime::g_lump_state[s];
        obs.push_back(L.x);
        obs.push_back(L.y);
        obs.push_back(L.z);
        obs.push_back(L.size);
        obs.push_back(L.mass);
        obs.push_back(L.peak);
        obs.push_back(L.min_lapse);
        obs.push_back(L.min_chi);
    }
    return obs;
}

inline double compute_l2_hamiltonian_norm(GRAMR &amr)
{
    auto &level0 = dynamic_cast<RadialRecipeLevel &>(amr.getLevel(0));
    return compute_radial_recipe_constraint_norms(level0).L2_Ham;
}

#endif /* RL_OBSERVATION_COLLECTOR_HPP_ */
