/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef AHFINDERLEVEL_HPP_
#define AHFINDERLEVEL_HPP_

#include "AlgebraicConstraintsEnforcer.hpp"
#include "BinaryBHInitialData.hpp"
#include "CCZ4RHS.hpp"
#include "DefaultLevelBld.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRAmr.hpp"
#include "GRAmrLevel.hpp"
#include "MovingPunctureGauge.hpp"
#include "PositiveChiAndLapse.hpp"

// We basically need this to have a valid AMR hierarchy. Initial data is the
// BinaryBH example's puncture data, evolved with CCZ4 and the moving puncture
// gauge as in the BinaryBH example. Tagging is a no-op.

class AHFinderLevel : public GRAmrLevel
{
  public:
    friend class DefaultLevelBld<AHFinderLevel>;

    // Inherit the contructors from GRAmrLevel
    using GRAmrLevel::GRAmrLevel;

    static void variableSetUp() { state_variable_set_up(); }

    // BinaryBH puncture initial data
    void initData() override
    {
        amrex::Real dx = Geom().CellSize(0);
        BinaryBHInitialData binary_initial_data(dx);

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

    // Enforce det(h)=1, trace free A_ij and positive chi and lapse after
    // each step
    void specific_advance() override
    {
        amrex::MultiFab &state_new = get_new_data(state_index);
        const auto &state_arrays   = state_new.arrays();

        AlgebraicConstraintsEnforcer algebraic_constraints_enforcer;
        PositiveChiAndLapse positive_chi_lapse;

        amrex::ParallelFor(
            state_new,
            [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
            {
                algebraic_constraints_enforcer(ix, iy, iz,
                                               state_arrays[box_no]);
                positive_chi_lapse(ix, iy, iz, state_arrays[box_no]);
            });
    }

    // CCZ4 with the moving puncture gauge, as in BinaryBHLevel
    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
    void specific_eval_rhs(amrex::MultiFab &a_soln, amrex::MultiFab &a_rhs,
                           const amrex::Real /*a_time*/) override
    {
        const auto &soln_arrays       = a_soln.arrays();
        const auto &const_soln_arrays = a_soln.const_arrays();
        const auto &rhs_arrays        = a_rhs.arrays();
        const auto soln_ghosts        = a_soln.nGrowVect();

        AlgebraicConstraintsEnforcer algebraic_constraints_enforcer;
        PositiveChiAndLapse positive_chi_lapse;

        amrex::ParallelFor(
            a_soln, soln_ghosts,
            [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
            {
                algebraic_constraints_enforcer(ix, iy, iz, soln_arrays[box_no]);
                positive_chi_lapse(ix, iy, iz, soln_arrays[box_no]);
            });

        CCZ4RHS<FourthOrderDerivatives> ccz4rhs(Geom().CellSize(0));
        MovingPunctureGauge<FourthOrderDerivatives> moving_puncture_gauge(
            Geom().CellSize(0));

        amrex::ParallelFor(
            a_rhs,
            [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
            {
                ccz4rhs.compute_chi_and_h_ij(ix, iy, iz, rhs_arrays[box_no],
                                             const_soln_arrays[box_no]);
            });

        amrex::ParallelFor(
            a_rhs,
            [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
            {
                ccz4rhs.compute_A_ij_and_Theta_and_Gamma(
                    ix, iy, iz, rhs_arrays[box_no], const_soln_arrays[box_no]);
            });

        amrex::ParallelFor(
            a_rhs,
            [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
            {
                moving_puncture_gauge.calculate_rhs(
                    ix, iy, iz, rhs_arrays[box_no], const_soln_arrays[box_no]);

                ccz4rhs.apply_dissipation(ix, iy, iz, rhs_arrays[box_no],
                                          const_soln_arrays[box_no]);
            });

        amrex::Gpu::streamSynchronize();
    }

    // Enforce det(h)=1 and trace free A_ij during the RK4 substeps
    void specific_update_ode(amrex::MultiFab &a_soln) override
    {
        AlgebraicConstraintsEnforcer algebraic_constraints_enforcer;

        const auto &soln_arrays = a_soln.arrays();
        amrex::ParallelFor(
            a_soln, amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
            {
                algebraic_constraints_enforcer(ix, iy, iz, soln_arrays[box_no]);
            });

        amrex::Gpu::streamSynchronize();
    }

    // No refinement needed
    void tag_cells(amrex::TagBoxArray &a_tag_box_array,
                   amrex::Real a_regrid_threshold) final
    {
    }
};

#endif /* AHFINDERLEVEL_HPP_ */
