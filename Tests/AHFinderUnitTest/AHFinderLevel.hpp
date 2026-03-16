/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef AHFINDERLEVEL_HPP_
#define AHFINDERLEVEL_HPP_

#include "DefaultLevelFactory.hpp"
#include "FixedGridsTagger.hpp"
#include "GRAMR.hpp"
#include "GRAMRLevel.hpp"

// We basically need this to have a valid AMR hierarchy

class AHFinderLevel : public GRAMRLevel
{
  public:
    friend class DefaultLevelFactory<AHFinderLevel>;

    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    static void variableSetUp() { stateVariableSetUp(); }

    // initialize data
    void initData()
    {
        amrex::MultiFab &state = get_new_data(State_Type);
        const auto &arrs       = state.arrays();
        auto const &geom       = Geom();
        auto const prob_lo     = geom.ProbLoArray();
        auto const dx          = geom.CellSizeArray();
        const int c_polystate  = 0; // index

        std::array<double, AMREX_SPACEDIM> center{AMREX_D_DECL(0., 0., 0.)};
        GRParmParse pp;
        pp.query("center", center);

        // Fill the state
        amrex::ParallelFor(
            state, state.nGrowVect(),
            [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
            {
                const auto &array = arrs[box_no];

                // compute coordinates
                amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0] - center[0];
                amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1] - center[1];
                amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2] - center[2];

                amrex::Real centre_distance = sqrt(x * x + y * y + z * z);

                // zero out everything first
                array(i, j, k, c_polystate) = 0.0;

                // write in
                array(i, j, k, c_polystate) = pow(centre_distance, 3) - 8;
            });

        amrex::Gpu::streamSynchronize();
    }

    void specificEvalRHS(amrex::MultiFab &a_soln, amrex::MultiFab &a_rhs,
                         const double a_time)
    {
    }

    void tag_cells(amrex::TagBoxArray &a_tag_box_array,
                   amrex::Real a_regrid_threshold) final
    {
        amrex::MultiFab &state_new = get_new_data(State_Type);

        const auto &tag_arrs   = a_tag_box_array.arrays();
        const auto &state_arrs = state_new.arrays();

        const amrex::Real dx         = Geom().CellSize(0);
        const int current_level      = Level();
        const amrex::Real box_length = Geom().ProbLength(0);
        std::array<double, AMREX_SPACEDIM> center{AMREX_D_DECL(0., 0., 0.)};
        GRParmParse pp;
        pp.query("center", center);

        FixedGridsTagger my_tagging_criterion{dx, current_level, box_length,
                                              center};

        amrex::ParallelFor(
            a_tag_box_array,
            [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
            { my_tagging_criterion(i, j, k, tag_arrs[box_no]); });
        amrex::Gpu::streamSynchronize();
    }
};

#endif /* AHFINDERLEVEL_HPP_ */