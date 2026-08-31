/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef PARTICLEINTERPOLATORLEVEL_HPP_
#define PARTICLEINTERPOLATORLEVEL_HPP_

#include "DefaultLevelBld.hpp"
#include "FixedGridsTagger.hpp"
#include "GRAmr.hpp"
#include "GRAmrLevel.hpp"
#include "PolynomialDerivedQuantity.hpp"
#include "StateTypes.hpp"

// We basically need this to have a valid AMR hierarchy

class ParticleInterpolatorLevel : public GRAmrLevel
{
  public:
    friend class DefaultLevelBld<ParticleInterpolatorLevel>;

    // Inherit the contructors from GRAmrLevel
    using GRAmrLevel::GRAmrLevel;

    static void variableSetUp()
    {
        state_variable_set_up();
        PolynomialDerivedQuantity::set_up(state_index);
    }

    // initialize data
    void initData() override
    {
        amrex::MultiFab &state = get_new_data(state_index);
        const auto &arrs       = state.arrays();
        auto const &geom       = Geom();
        auto const prob_lo     = geom.ProbLoArray();
        auto const dx          = geom.CellSizeArray();

        std::array<amrex::Real, AMREX_SPACEDIM> center{
            AMREX_D_DECL(0., 0., 0.)};
        GRParmParse pp;
        pp.query("geometry.center", center);

        // Fill the state
        amrex::ParallelFor(
            state, state.nGrowVect(),
            // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
            [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
            {
                const auto &array = arrs[box_no];

                // compute coordinates
                amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2] - center[2];

                // write in
                array(i, j, k, c_polystate) = z * z * z;
            });

        amrex::Gpu::streamSynchronize();
    }

    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
    void specific_eval_rhs(amrex::MultiFab &a_soln, amrex::MultiFab &a_rhs,
                           const amrex::Real a_time) override
    {
    }

    void tag_cells(amrex::TagBoxArray &a_tag_box_array,
                   amrex::Real a_regrid_threshold) final
    {
        amrex::MultiFab &state_new = get_new_data(state_index);

        const auto &tag_arrs   = a_tag_box_array.arrays();
        const auto &state_arrs = state_new.arrays();

        const amrex::Real dx    = Geom().CellSize(0);
        const int current_level = Level();

        FixedGridsTagger my_tagging_criterion{dx, current_level};

        // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
        amrex::ParallelFor(
            a_tag_box_array,
            [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
            { my_tagging_criterion(i, j, k, tag_arrs[box_no]); });
        amrex::Gpu::streamSynchronize();
    }
};

#endif /* PARTICLEINTERPOLATORLEVEL_HPP_ */
