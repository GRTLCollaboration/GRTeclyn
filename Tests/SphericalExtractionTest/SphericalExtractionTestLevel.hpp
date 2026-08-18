/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SPHERICALEXTRACTIONTESTLEVEL_HPP_
#define SPHERICALEXTRACTIONTESTLEVEL_HPP_

#include "DefaultLevelFactory.hpp"
#include "FixedGridsTagger.hpp"
#include "GRAMR.hpp"
#include "GRAMRLevel.hpp"
#include "GRParmParse.hpp"
#include "SimulationParameters.hpp"
#include "SphericalHarmonics.hpp"
#include "StateTypes.hpp"
#include "StateVariables.hpp"

class SphericalExtractionTestLevel : public GRAMRLevel
{
  public:
    friend class DefaultLevelFactory<SphericalExtractionTestLevel>;

    using GRAMRLevel::GRAMRLevel;

    static void variableSetUp() { stateVariableSetUp(); }

    // save data in state vars
    void initData() override
    {
        amrex::MultiFab &state = get_new_data(state_index);
        const auto &arrs       = state.arrays();
        auto const &geom       = Geom();
        auto const prob_lo     = geom.ProbLoArray();
        auto const dx          = geom.CellSizeArray();

        GRParmParse pp;
        GRParmParse test_pp("test");

        int es{0};
        int el{2};
        int em{0}; // spherical harmonic params
        test_pp.get("es", es);
        test_pp.get("el", el);
        test_pp.get("em", em);

        std::array<double, AMREX_SPACEDIM> center{};
        pp.get("geometry.center", center);

        // Fill the state
        amrex::ParallelFor(
            state, state.nGrowVect(),
            // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
            [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
            {
                const auto &array = arrs[box_no];

                const amrex::Real x =
                    prob_lo[0] + (i + 0.5) * dx[0] - center[0];
                const amrex::Real y =
                    prob_lo[1] + (j + 0.5) * dx[1] - center[1];
                const amrex::Real z =
                    prob_lo[2] + (k + 0.5) * dx[2] - center[2];

                const auto Y_lm =
                    SphericalHarmonics::spin_Y_lm(x, y, z, es, el, em);

                array(i, j, k, c_phi_Re) = Y_lm.Real;
                array(i, j, k, c_phi_Im) = Y_lm.Im;
            });

        amrex::Gpu::streamSynchronize();
    }

    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
    void specificEvalRHS(amrex::MultiFab &a_soln, amrex::MultiFab &a_rhs,
                         const double a_time) override
    {
    }

    void tag_cells(amrex::TagBoxArray &a_tag_box_array,
                   amrex::Real a_regrid_threshold) final
    {
        amrex::MultiFab &state_new = get_new_data(state_index);

        const auto &tag_arrs   = a_tag_box_array.arrays();
        const auto &state_arrs = state_new.arrays();

        const amrex::Real dx         = Geom().CellSize(0);
        const int current_level      = Level();
        const amrex::Real box_length = Geom().ProbLength(0);

        std::array<double, AMREX_SPACEDIM> center{AMREX_D_DECL(0., 0., 0.)};
        GRParmParse pp;
        pp.query("geometry.center", center);

        FixedGridsTagger my_tagging_criterion{dx, current_level, box_length,
                                              center};

        // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
        amrex::ParallelFor(
            a_tag_box_array,
            [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
            { my_tagging_criterion(i, j, k, tag_arrs[box_no]); });

        amrex::Gpu::streamSynchronize();
    }
};

#endif /* SPHERICALEXTRACTIONTESTLEVEL_HPP_ */
