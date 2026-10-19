/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#include "ScalarFieldLevel.hpp"

#include "AMRReductions.hpp"
#include "AlgebraicConstraintsEnforcer.hpp"
#include "CCZ4RHSWithMatter.hpp"
#include "ConstraintsWithMatter.hpp"
#include "EMTensor.hpp"
#include "FilesystemTools.hpp"
#include "FixedGridsTagger.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRParmParse.hpp"
#include "GammaCalculator.hpp"
#include "IntegratedMovingPunctureGauge.hpp"
#include "LineExtraction.hpp"
#include "OscillatonInitialData.hpp"
#include "PositiveChiAndLapse.hpp"
#include "SixthOrderDerivatives.hpp"
#include "SmallDataIO.hpp"
#include "StateTypes.hpp"

#include <type_traits>

using ScalarFieldEnergyDensity =
    EMTensor<ScalarFieldLevel::ScalarFieldWithPotential<>,
             EMTensorOptions::justEnergyDensity>;
using ScalarFieldConstraints =
    ConstraintsWithMatter<ScalarFieldLevel::ScalarFieldWithPotential<>>;

ScalarFieldAmr *ScalarFieldLevel::get_scalar_field_amr_ptr()
{
    return dynamic_cast<ScalarFieldAmr *>(get_gr_amr_ptr());
}

void ScalarFieldLevel::variableSetUp()
{
    BL_PROFILE("ScalarFieldLevel::variableSetUp()");
    state_variable_set_up();
    ScalarFieldConstraints::set_up(state_index);
    ScalarFieldEnergyDensity::set_up(state_index);
}

void ScalarFieldLevel::specific_advance()
{
    BL_PROFILE("ScalarFieldLevel::specific_advance()");

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

void ScalarFieldLevel::specific_post_timestep()
{
    BL_PROFILE("ScalarFieldLevel::specific_post_timestep()");

    if (Level() == 0)
    {
        const amrex::Real time         = get_state_data(state_index).curTime();
        const amrex::Real dt           = get_gr_amr_ptr()->dtLevel(0);
        const amrex::Real restart_time = get_gr_amr_ptr()->get_restart_time();
        const bool first_step          = (time <= dt);

        const LineExtraction<1> phi_extraction("phi_line_extraction", c_phi, dt,
                                               time, restart_time, first_step);
        phi_extraction.execute_query(
            &get_scalar_field_amr_ptr()->phi_interpolator);

        const LineExtraction<1> rho_extraction("rho_line_extraction", 0, dt,
                                               time, restart_time, first_step);
        const LineExtraction<1>::derived_vars_t rho_vars{
            .name            = ScalarFieldEnergyDensity::name,
            .component_names = {"rho"},
            .parities        = {BCParity::even},
        };

        rho_extraction.execute_query(
            &get_scalar_field_amr_ptr()->rho_interpolator, rho_vars);

        GRParmParse pp;
        bool store_reductions = false;
        pp.query("diagnostic.store_reductions", store_reductions);
        if (store_reductions)
        {
            write_reductions(time, dt, restart_time, first_step);
        }
    }
}

void ScalarFieldLevel::write_reductions(amrex::Real a_time, amrex::Real a_dt,
                                        amrex::Real a_restart_time,
                                        bool a_first_step)
{
    BL_PROFILE("ScalarFieldLevel::write_reductions()");

    // Derive the constraints (Ham, Mom1-3) and phi into a single hierarchy,
    // so that all reductions are computed in one pass over the whole hierarchy.
    auto &amr = *get_gr_amr_ptr();
    const AMRReductions::DerivedHierarchy diagnostics(
        amr, a_time, {Constraints::name, "phi"});

    const AMRReductions::AMRReducer reducer(amr);
    const auto [ham_l2_norm, mom_l2_norm, phi_max_abs] =
        reducer.compute(diagnostics.data(),
                        AMRReductions::L2Norm{.comp = diagnostics.comp("Ham")},
                        AMRReductions::L2Norm{
                            .comp      = diagnostics.comp("Mom1"),
                            .num_comps = AMREX_SPACEDIM,
                        },
                        AMRReductions::MaxAbs{.comp = diagnostics.comp("phi")});

    GRParmParse pp;
    std::string output_path;
    pp.get("grteclyn.output_path", output_path);
    const std::string diagnostics_path = output_path + "/diagnostics";
    FilesystemTools::ensure_directory_exists(diagnostics_path);

    SmallDataIO constraints_file(diagnostics_path + "/constraint_norms", a_dt,
                                 a_time, a_restart_time, SmallDataIO::APPEND,
                                 a_first_step);
    SmallDataIO phi_max_abs_file(diagnostics_path + "/max_abs_phi", a_dt,
                                 a_time, a_restart_time, SmallDataIO::APPEND,
                                 a_first_step);
    if (a_first_step)
    {
        constraints_file.write_header_line({"L2_Ham", "L2_Mom"});
        phi_max_abs_file.write_header_line({"max_abs_phi"});
    }
    constraints_file.write_time_data_line(
        std::vector<amrex::Real>{ham_l2_norm, mom_l2_norm});
    phi_max_abs_file.write_time_data_line(
        std::vector<amrex::Real>{phi_max_abs});
}

void ScalarFieldLevel::initData()
{
    BL_PROFILE("ScalarFieldLevel::initData()");

    if (get_gr_amr_ptr()->Verbose() > 0)
    {
        amrex::Print() << "ScalarFieldLevel::initData " << Level() << "\n";
    }

    amrex::MultiFab &state_new = get_new_data(state_index);
    const auto &state_arrays   = state_new.arrays();

    const OscillatonInitialData initial_data(Geom().CellSize(0));
    static_assert(std::is_trivially_copyable_v<OscillatonInitialData>,
                  "OscillatonInitialData must be device copyable");

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
void ScalarFieldLevel::specific_eval_rhs(amrex::MultiFab &a_soln,
                                         amrex::MultiFab &a_rhs,
                                         const amrex::Real /*a_time*/)
{
    BL_PROFILE("ScalarFieldLevel::specific_eval_rhs()");

    const auto &soln_arrays       = a_soln.arrays();
    const auto &const_soln_arrays = a_soln.const_arrays();
    const auto &rhs_arrays        = a_rhs.arrays();

    const AlgebraicConstraintsEnforcer algebraic_constraints_enforcer;
    const PositiveChiAndLapse positive_chi_and_lapse;

    amrex::ParallelFor(
        a_soln, a_soln.nGrowVect(),
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
        const CCZ4RHSWithMatter<
            ScalarFieldWithPotential<FourthOrderDerivatives>,
            FourthOrderDerivatives>
            ccz4_rhs(Geom().CellSize(0));
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
                ccz4_rhs.add_emtensor_rhs(ix, iy, iz, rhs_arrays[box_no],
                                          const_soln_arrays[box_no]);
                integrated_moving_puncture_gauge.calculate_rhs(
                    ix, iy, iz, rhs_arrays[box_no], const_soln_arrays[box_no]);
                ccz4_rhs.add_matter_rhs(ix, iy, iz, rhs_arrays[box_no],
                                        const_soln_arrays[box_no]);
                ccz4_rhs.apply_dissipation(ix, iy, iz, rhs_arrays[box_no],
                                           const_soln_arrays[box_no]);
            });
    }
    else if (m_evolution_spatial_derivative_order == 6)
    {
        const CCZ4RHSWithMatter<ScalarFieldWithPotential<SixthOrderDerivatives>,
                                SixthOrderDerivatives>
            ccz4_rhs(Geom().CellSize(0));
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
                ccz4_rhs.add_emtensor_rhs(ix, iy, iz, rhs_arrays[box_no],
                                          const_soln_arrays[box_no]);
                integrated_moving_puncture_gauge.calculate_rhs(
                    ix, iy, iz, rhs_arrays[box_no], const_soln_arrays[box_no]);
                ccz4_rhs.add_matter_rhs(ix, iy, iz, rhs_arrays[box_no],
                                        const_soln_arrays[box_no]);
                ccz4_rhs.apply_dissipation(ix, iy, iz, rhs_arrays[box_no],
                                           const_soln_arrays[box_no]);
            });
    }
    // NOLINTEND(bugprone-branch-clone)

    amrex::Gpu::streamSynchronize();
}

void ScalarFieldLevel::specific_update_ode(amrex::MultiFab &a_soln)
{
    BL_PROFILE("ScalarFieldLevel::specific_update_ode()");

    const auto &soln_arrays = a_soln.arrays();
    const AlgebraicConstraintsEnforcer algebraic_constraints_enforcer;

    amrex::ParallelFor(
        a_soln, amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
        { algebraic_constraints_enforcer(ix, iy, iz, soln_arrays[box_no]); });
    amrex::Gpu::streamSynchronize();
}

void ScalarFieldLevel::tag_cells(amrex::TagBoxArray &a_tag_box_array,
                                 const amrex::Real /*a_regrid_threshold*/)
{
    BL_PROFILE("ScalarFieldLevel::tag_cells()");

    const auto &tag_arrays = a_tag_box_array.arrays();

    const FixedGridsTagger tagger(Geom().CellSize(0), Level());

    amrex::ParallelFor(a_tag_box_array,
                       [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
                       { tagger(ix, iy, iz, tag_arrays[box_no]); });
    amrex::Gpu::streamSynchronize();
}
