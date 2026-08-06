/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef AHFINDERLEVEL_HPP_
#define AHFINDERLEVEL_HPP_

#include "BinaryBHInitialData.hpp"
#include "DefaultLevelFactory.hpp"
#include "GRAMR.hpp"
#include "GRAMRLevel.hpp"

// We basically need this to have a valid AMR hierarchy. Initial data is the
// BinaryBH example's puncture data -- we only need the first timestep (no
// evolution), so this level does not implement specificEvalRHS or tagging
// beyond a no-op.

class AHFinderLevel : public GRAMRLevel
{
  public:
    friend class DefaultLevelFactory<AHFinderLevel>;

    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    static void variableSetUp() { stateVariableSetUp(); }

    // BinaryBH puncture initial data
    void initData()
    {
        double dx = Geom().CellSize(0);
        BinaryBHInitialData binary_initial_data(simParams().bh1_params,
                                                simParams().bh2_params, dx);

        static_assert(std::is_trivially_copyable_v<BinaryBHInitialData>,
                      "BinaryBHInitialData needs to be device copyable");
        amrex::MultiFab &state_new = get_new_data(state_index);
        const auto &state_arrays   = state_new.arrays();
        amrex::ParallelFor(
            state_new, state_new.nGrowVect(),
            [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
            {
                amrex::CellData<amrex::Real> cell =
                    state_arrays[box_no].cellData(ix, iy, iz);
                for (int n = 0; n < cell.nComp(); ++n)
                {
                    cell[n] = 0.;
                }
                binary_initial_data(ix, iy, iz, state_arrays[box_no]);
            });

        amrex::Gpu::streamSynchronize();
    }

    void specificEvalRHS(amrex::MultiFab &a_soln, amrex::MultiFab &a_rhs,
                         const double a_time)
    {
    }

    // No refinement needed
    void tag_cells(amrex::TagBoxArray &a_tag_box_array,
                   amrex::Real a_regrid_threshold) final
    {
    }
};

#endif /* AHFINDERLEVEL_HPP_ */
