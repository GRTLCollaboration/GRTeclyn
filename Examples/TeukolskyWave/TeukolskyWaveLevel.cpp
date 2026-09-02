/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#include "TeukolskyWaveLevel.hpp"

#include "AlgebraicConstraintsEnforcer.hpp"
#include "CCZ4RHS.hpp"
#include "Constraints.hpp"
#include "EppleyPacket.hpp"
#include "FixedGridsTagger.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GammaCalculator.hpp"
#include "IntegratedMovingPunctureGauge.hpp"
#include "PositiveChiAndLapse.hpp"
#include "SixthOrderDerivatives.hpp"
#include "StateTypes.hpp"
#include "TeukolskyWaveInitialData.hpp"
#include "Weyl4.hpp"
#include "WeylExtraction.hpp"

void TeukolskyWaveLevel::variableSetUp()
{
    BL_PROFILE("TeukolskyWaveLevel::variableSetUp()");
    state_variable_set_up();
    Constraints::set_up(state_index);
    Weyl4::set_up(state_index);
}

void TeukolskyWaveLevel::specific_advance()
{
    BL_PROFILE("TeukolskyWaveLevel::specific_advance()");

    amrex::MultiFab &state_new = get_new_data(state_index);
    const auto &state_arrays   = state_new.arrays();

    const AlgebraicConstraintsEnforcer algebraic_constraints_enforcer;
    const PositiveChiAndLapse positive_chi_and_lapse;

    amrex::ParallelFor(
        state_new,
        [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
        {
            algebraic_constraints_enforcer(ix, iy, iz, state_arrays[box_no]);
            positive_chi_and_lapse(ix, iy, iz, state_arrays[box_no]);
        });
    amrex::Gpu::streamSynchronize();
}

void TeukolskyWaveLevel::specific_post_timestep()
{
    BL_PROFILE("TeukolskyWaveLevel::specific_post_timestep()");

    spherical_extraction_params_t extraction_params("weyl_extraction");
    extraction_params.fill_params();

    if (extraction_params.enabled)
    {
        const int min_level = extraction_params.min_extraction_level();
        bool calculate_weyl = at_level_timestep_multiple(min_level);

        if (calculate_weyl && Level() == min_level)
        {
            amrex::Real m_time       = get_state_data(state_index).curTime();
            amrex::Real m_dt         = get_gr_amr_ptr()->dtLevel(Level());
            amrex::Real restart_time = get_gr_amr_ptr()->get_restart_time();
            bool first_step          = (m_time <= m_dt);

            WeylExtraction my_extraction(extraction_params, m_dt, m_time,
                                         first_step, restart_time);
            my_extraction.execute_query(&get_gr_amr_ptr()->m_weyl_interpolator);
        }
    }
}

void TeukolskyWaveLevel::initData()
{
    BL_PROFILE("TeukolskyWaveLevel::initData()");

    if (get_gr_amr_ptr()->Verbose() > 0)
    {
        amrex::Print() << "TeukolskyWaveLevel::initData " << Level() << "\n";
    }

    // Load the parity and magnetic number from the parameters file.
    int magnetic{};
    std::string parity{};
    amrex::ParmParse pp;
    pp.get("teukolsky_wave.magnetic", magnetic);
    pp.get("teukolsky_wave.parity", parity);

    amrex::MultiFab &state_new = get_new_data(state_index);
    const auto &state_arrays   = state_new.arrays();

    const TeukolskyWaveInitialData initial_data(Geom().CellSize(0));
    initial_data.initialize_eppley_packet(magnetic, parity);

    static_assert(std::is_trivially_copyable_v<TeukolskyWaveInitialData>,
                  "TeukolskyWaveInitialData must be device copyable");

    amrex::ParallelFor(state_new, state_new.nGrowVect(),
                       [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
                       {
                           const amrex::CellData<amrex::Real> cell =
                               state_arrays[box_no].cellData(ix, iy, iz);
                           for (int component = 0; component < cell.nComp();
                                ++component)
                           {
                               cell[component] = 0.0;
                           }
                           initial_data(ix, iy, iz, state_arrays[box_no]);
                       });

    if (m_evolution_spatial_derivative_order == 4)
    {
        const GammaCalculator<FourthOrderDerivatives> gamma_calculator(
            Geom().CellSize(0));
        const IntegratedMovingPunctureGauge<FourthOrderDerivatives>
            integrated_moving_puncture_gauge(Geom().CellSize(0));
        amrex::ParallelFor(
            state_new,
            [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
            {
                gamma_calculator(ix, iy, iz, state_arrays[box_no]);
                integrated_moving_puncture_gauge.set_initial_B_to_Gamma(
                    ix, iy, iz, state_arrays[box_no]);
            });
    }
    else if (m_evolution_spatial_derivative_order == 6)
    {
        const GammaCalculator<SixthOrderDerivatives> gamma_calculator(
            Geom().CellSize(0));
        const IntegratedMovingPunctureGauge<SixthOrderDerivatives>
            integrated_moving_puncture_gauge(Geom().CellSize(0));
        amrex::ParallelFor(
            state_new,
            [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
            {
                gamma_calculator(ix, iy, iz, state_arrays[box_no]);
                integrated_moving_puncture_gauge.set_initial_B_to_Gamma(
                    ix, iy, iz, state_arrays[box_no]);
            });
    }

    amrex::Gpu::streamSynchronize();
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
void TeukolskyWaveLevel::specific_eval_rhs(amrex::MultiFab &a_soln,
                                           amrex::MultiFab &a_rhs,
                                           const amrex::Real /*a_time*/)
{
    BL_PROFILE("TeukolskyWaveLevel::specific_eval_rhs()");

    const auto &soln_arrays       = a_soln.arrays();
    const auto &const_soln_arrays = a_soln.const_arrays();
    const auto &rhs_arrays        = a_rhs.arrays();
    const auto soln_ghosts        = a_soln.nGrowVect();

    const AlgebraicConstraintsEnforcer algebraic_constraints_enforcer;
    const PositiveChiAndLapse positive_chi_and_lapse;

    amrex::ParallelFor(
        a_soln, soln_ghosts,
        [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
        {
            algebraic_constraints_enforcer(ix, iy, iz, soln_arrays[box_no]);
            positive_chi_and_lapse(ix, iy, iz, soln_arrays[box_no]);
        });

    if (m_evolution_spatial_derivative_order != 4 &&
        m_evolution_spatial_derivative_order != 6)
    {
        amrex::Abort("spatial_derivative_order must be 4 or 6");
    }

    // NOLINTBEGIN(bugprone-branch-clone)
    if (m_evolution_spatial_derivative_order == 4)
    {
        const CCZ4RHS<FourthOrderDerivatives> ccz4_rhs(Geom().CellSize(0));
        const IntegratedMovingPunctureGauge<FourthOrderDerivatives>
            integrated_moving_puncture_gauge(Geom().CellSize(0));

        amrex::ParallelFor(
            a_rhs,
            [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
            {
                ccz4_rhs.compute_chi_and_h_ij(ix, iy, iz, rhs_arrays[box_no],
                                              const_soln_arrays[box_no]);
            });

        // NOLINTBEGIN(bugprone-easily-swappable-parameters)
        amrex::ParallelFor(
            a_rhs,
            [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
            {
                ccz4_rhs.compute_A_ij_and_Theta_and_Gamma(
                    ix, iy, iz, rhs_arrays[box_no], const_soln_arrays[box_no]);
            });
        // NOLINTEND(bugprone-easily-swappable-parameters)

        amrex::ParallelFor(
            a_rhs,
            [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
            {
                integrated_moving_puncture_gauge.calculate_rhs(
                    ix, iy, iz, rhs_arrays[box_no], const_soln_arrays[box_no]);
                ccz4_rhs.apply_dissipation(ix, iy, iz, rhs_arrays[box_no],
                                           const_soln_arrays[box_no]);
            });
    }
    else if (m_evolution_spatial_derivative_order == 6)
    {
        const CCZ4RHS<SixthOrderDerivatives> ccz4_rhs(Geom().CellSize(0));
        const IntegratedMovingPunctureGauge<SixthOrderDerivatives>
            integrated_moving_puncture_gauge(Geom().CellSize(0));

        amrex::ParallelFor(
            a_rhs,
            [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
            {
                ccz4_rhs.compute_chi_and_h_ij(ix, iy, iz, rhs_arrays[box_no],
                                              const_soln_arrays[box_no]);
            });

        // NOLINTBEGIN(bugprone-easily-swappable-parameters)
        amrex::ParallelFor(
            a_rhs,
            [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
            {
                ccz4_rhs.compute_A_ij_and_Theta_and_Gamma(
                    ix, iy, iz, rhs_arrays[box_no], const_soln_arrays[box_no]);
            });
        // NOLINTEND(bugprone-easily-swappable-parameters)

        amrex::ParallelFor(
            a_rhs,
            [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
            {
                integrated_moving_puncture_gauge.calculate_rhs(
                    ix, iy, iz, rhs_arrays[box_no], const_soln_arrays[box_no]);
                ccz4_rhs.apply_dissipation(ix, iy, iz, rhs_arrays[box_no],
                                           const_soln_arrays[box_no]);
            });
    }
    // NOLINTEND(bugprone-branch-clone)

    amrex::Gpu::streamSynchronize();
}

void TeukolskyWaveLevel::specific_update_ode(amrex::MultiFab &a_soln)
{
    BL_PROFILE("TeukolskyWaveLevel::specific_update_ode()");

    const auto &soln_arrays = a_soln.arrays();
    const AlgebraicConstraintsEnforcer algebraic_constraints_enforcer;

    amrex::ParallelFor(
        a_soln, amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
        { algebraic_constraints_enforcer(ix, iy, iz, soln_arrays[box_no]); });
    amrex::Gpu::streamSynchronize();
}

void TeukolskyWaveLevel::tag_cells(amrex::TagBoxArray &a_tag_box_array,
                                   const amrex::Real /*a_regrid_threshold*/)
{
    BL_PROFILE("TeukolskyWaveLevel::tag_cells()");

    const auto &tag_arrays = a_tag_box_array.arrays();

    const FixedGridsTagger tagger(Geom().CellSize(0), Level());

    amrex::ParallelFor(a_tag_box_array,
                       [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
                       { tagger(ix, iy, iz, tag_arrays[box_no]); });
    amrex::Gpu::streamSynchronize();
}
