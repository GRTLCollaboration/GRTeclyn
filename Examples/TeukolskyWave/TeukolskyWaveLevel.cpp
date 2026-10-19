/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#include "TeukolskyWaveLevel.hpp"

#include "AMRReductions.hpp"
#include "AlgebraicConstraintsEnforcer.hpp"
#include "CCZ4RHS.hpp"
#include "Constraints.hpp"
#include "EppleyPacket.hpp"
#include "FilesystemTools.hpp"
#include "FixedGridsTagger.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRParmParse.hpp"
#include "GammaCalculator.hpp"
#include "IntegratedMovingPunctureGauge.hpp"
#include "PositiveChiAndLapse.hpp"
#include "SixthOrderDerivatives.hpp"
#include "SmallDataIO.hpp"
#include "StateTypes.hpp"
#include "TeukolskyWaveInitialData.hpp"
#include "Weyl4.hpp"
#include "WeylExtraction.hpp"

TeukolskyWaveAmr *TeukolskyWaveLevel::get_tw_amr_ptr()
{
    return dynamic_cast<TeukolskyWaveAmr *>(get_gr_amr_ptr());
}

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
            my_extraction.execute_query(&get_tw_amr_ptr()->m_weyl_interpolator);
        }
    }
    GRParmParse pp;
    bool store_constraint_norms = false;
    pp.query("diagnostic.store_constraint_norms", store_constraint_norms);
    if (store_constraint_norms && Level() == 0)
    {
        const amrex::Real time         = get_state_data(state_index).curTime();
        const amrex::Real dt           = get_gr_amr_ptr()->dtLevel(Level());
        const amrex::Real restart_time = get_gr_amr_ptr()->get_restart_time();
        const bool first_step          = (time <= dt);

        // Derive the constraints (Ham, Mom1-3) and chi into a single
        // hierarchy, so that all reductions below are computed in one pass
        // over the grid.
        auto &amr = *get_gr_amr_ptr();
        const AMRReductions::DerivedHierarchy diagnostics(
            amr, time, {Constraints::name, "chi"});

        // Optionally exclude a layer of cells next to the faces with a
        // Sommerfeld boundary condition, e.g. to check whether the cells whose
        // evolution equation it replaces dominate the norm. Reflective faces
        // are symmetry planes, evolved like the interior, so are kept.
        int boundary_cells_to_exclude = 0;
        pp.query("diagnostic.constraint_norms_boundary_cells_to_exclude",
                 boundary_cells_to_exclude);
        AMRReductions::BoundaryExclusion exclusion;
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
        {
            const auto is_sommerfeld = [&](amrex::Orientation::Side a_side)
            {
                return m_boundaries.get_boundary_condition(amrex::Orientation(
                           dir, a_side)) == BoundaryConditions::SOMMERFELD_BC;
            };
            if (is_sommerfeld(amrex::Orientation::low))
            {
                exclusion.lo[dir] = boundary_cells_to_exclude;
            }
            if (is_sommerfeld(amrex::Orientation::high))
            {
                exclusion.hi[dir] = boundary_cells_to_exclude;
            }
        }

        // The momentum constraint components Mom1-3 (which are consecutive)
        // are combined into a single Mom L2 norm.
        const AMRReductions::AMRReducer reducer(amr, exclusion);
        const auto [ham_l2_norm, mom_l2_norm, chi_min] = reducer.compute(
            diagnostics.data(),
            AMRReductions::L2Norm{.comp = diagnostics.comp("Ham")},
            AMRReductions::L2Norm{.comp      = diagnostics.comp("Mom1"),
                                  .num_comps = AMREX_SPACEDIM},
            AMRReductions::Min{.comp = diagnostics.comp("chi")});

        std::string output_path;
        pp.get("grteclyn.output_path", output_path);
        const std::string diagnostics_path = output_path + "/diagnostics";
        FilesystemTools::ensure_directory_exists(diagnostics_path);

        SmallDataIO constraints_file(diagnostics_path + "/constraint_norms", dt,
                                     time, restart_time, SmallDataIO::APPEND,
                                     first_step);
        SmallDataIO min_chi_file(diagnostics_path + "/min_chi", dt, time,
                                 restart_time, SmallDataIO::APPEND, first_step);
        if (first_step)
        {
            constraints_file.write_header_line({"L2_Ham", "L2_Mom"});
            min_chi_file.write_header_line({"min_chi"});
        }
        constraints_file.write_time_data_line(
            std::vector<amrex::Real>{ham_l2_norm, mom_l2_norm});
        min_chi_file.write_time_data_line(std::vector<amrex::Real>{chi_min});
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

    TeukolskyWaveInitialData initial_data(Geom().CellSize(0));
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
