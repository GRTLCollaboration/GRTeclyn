/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#include "ScalarFieldLevel.hpp"

#include "CCZ4RHSWithMatter.hpp"
#include "EMTensor.hpp"
#include "FixedGridsTagger.hpp"
#include "GammaCalculator.hpp"
#include "LineExtraction.hpp"
#include "MovingPunctureGaugeWithMatter.hpp"
#include "OscillatonInitialData.hpp"
#include "PositiveChiAndLapse.hpp"
#include "StateTypes.hpp"
#include "TraceARemoval.hpp"

#include <cmath>
#include <type_traits>

using ScalarFieldEnergyDensity =
    EMTensor<ScalarFieldLevel::ScalarFieldWithPotential,
             EMTensorOptions::justEnergyDensity>;

ScalarFieldAMR *ScalarFieldLevel::get_scalar_field_amr_ptr()
{
    return dynamic_cast<ScalarFieldAMR *>(get_gramr_ptr());
}

void ScalarFieldLevel::variableSetUp()
{
    BL_PROFILE("ScalarFieldLevel::variableSetUp()");
    stateVariableSetUp();
    ScalarFieldEnergyDensity::set_up(state_index);
}

void ScalarFieldLevel::specificAdvance()
{
    BL_PROFILE("ScalarFieldLevel::specificAdvance()");

    amrex::MultiFab &state_new = get_new_data(state_index);
    const auto &state_arrays   = state_new.arrays();

    const TraceARemoval trace_A_removal;
    const PositiveChiAndLapse positive_chi_and_lapse;

    amrex::ParallelFor(state_new,
                       [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
                       {
                           trace_A_removal(ix, iy, iz, state_arrays[box_no]);
                           positive_chi_and_lapse(ix, iy, iz,
                                                  state_arrays[box_no]);
                       });
    amrex::Gpu::streamSynchronize();
}

void ScalarFieldLevel::specificPostTimeStep()
{
    BL_PROFILE("ScalarFieldLevel::specificPostTimeStep()");

    if (Level() == 0 && simParams().activate_line_extraction)
    {
        const amrex::Real time         = get_state_data(state_index).curTime();
        const amrex::Real dt           = get_gramr_ptr()->dtLevel(0);
        const amrex::Real restart_time = get_gramr_ptr()->get_restart_time();
        const bool first_step          = (time <= dt);

        const auto start_coords = simParams().center;
        auto end_coords         = start_coords;
        const amrex::Real coordinate_offset =
            simParams().line_extraction_max_radius /
            std::sqrt(static_cast<amrex::Real>(AMREX_SPACEDIM));
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
        {
            end_coords[dir] += coordinate_offset;
        }

        const LineExtraction<1> phi_extraction(
            c_phi, simParams().line_extraction_num_points, start_coords,
            end_coords, dt, time, restart_time, first_step);
        phi_extraction.execute_query(
            &get_scalar_field_amr_ptr()->phi_interpolator,
            simParams().data_path + "phi_profile_");

        const LineExtraction<1> rho_extraction(
            0, simParams().line_extraction_num_points, start_coords, end_coords,
            dt, time, restart_time, first_step);
        const LineExtraction<1>::derived_vars_t rho_vars{
            ScalarFieldEnergyDensity::name, {"rho"}, {BCParity::even}};
        rho_extraction.execute_query(
            &get_scalar_field_amr_ptr()->rho_interpolator,
            simParams().data_path + "rho_profile_", rho_vars);
    }
}

void ScalarFieldLevel::initData()
{
    BL_PROFILE("ScalarFieldLevel::initData()");

    if (m_verbosity > 0)
    {
        amrex::Print() << "ScalarFieldLevel::initData " << Level() << "\n";
    }

    amrex::MultiFab &state_new = get_new_data(state_index);
    const auto &state_arrays   = state_new.arrays();

    const OscillatonInitialData initial_data(simParams().initial_params,
                                             Geom().CellSize(0));
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

    const GammaCalculator gamma_calculator(Geom().CellSize(0));
    amrex::ParallelFor(state_new,
                       [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
                       { gamma_calculator(ix, iy, iz, state_arrays[box_no]); });

    amrex::Gpu::streamSynchronize();
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
void ScalarFieldLevel::specificEvalRHS(amrex::MultiFab &a_soln,
                                       amrex::MultiFab &a_rhs,
                                       const double /*a_time*/)
{
    BL_PROFILE("ScalarFieldLevel::specificEvalRHS()");

    const auto &soln_arrays       = a_soln.arrays();
    const auto &const_soln_arrays = a_soln.const_arrays();
    const auto &rhs_arrays        = a_rhs.arrays();

    const TraceARemoval trace_A_removal;
    const PositiveChiAndLapse positive_chi_and_lapse;

    amrex::ParallelFor(a_soln, a_soln.nGrowVect(),
                       [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
                       {
                           trace_A_removal(ix, iy, iz, soln_arrays[box_no]);
                           positive_chi_and_lapse(ix, iy, iz,
                                                  soln_arrays[box_no]);
                       });

    if (simParams().max_spatial_derivative_order != 4)
    {
        amrex::Abort("ScalarField currently supports fourth-order spatial "
                     "derivatives only");
    }

    const Potential potential(simParams().potential_params);
    const ScalarFieldWithPotential scalar_field(potential,
                                                simParams().G_Newton);
    const CCZ4RHSWithMatter<ScalarFieldWithPotential,
                            MovingPunctureGaugeWithMatter,
                            FourthOrderDerivatives>
        ccz4_rhs(scalar_field, simParams().ccz4_params, Geom().CellSize(0),
                 simParams().sigma, CCZ4RHS<>::USE_CCZ4);

    amrex::ParallelFor(a_rhs,
                       [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
                       {
                           ccz4_rhs.compute_chi_and_h_ij(
                               ix, iy, iz, rhs_arrays[box_no],
                               const_soln_arrays[box_no]);
                       });

    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    amrex::ParallelFor(
        a_rhs,
        [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
        {
            ccz4_rhs
                .compute_A_ij_and_Theta_and_Gamma<CCZ4RHS<>::USE_CCZ4, true>(
                    ix, iy, iz, rhs_arrays[box_no], const_soln_arrays[box_no]);
        });
    // NOLINTEND(bugprone-easily-swappable-parameters)

    amrex::ParallelFor(
        a_rhs,
        [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
        {
            ccz4_rhs.calculate_gauge_rhs(ix, iy, iz, rhs_arrays[box_no],
                                         const_soln_arrays[box_no]);
            ccz4_rhs.add_emtensor_rhs(ix, iy, iz, rhs_arrays[box_no],
                                      const_soln_arrays[box_no]);
            ccz4_rhs.add_matter_rhs(ix, iy, iz, rhs_arrays[box_no],
                                    const_soln_arrays[box_no]);
            ccz4_rhs.apply_dissipation(ix, iy, iz, rhs_arrays[box_no],
                                       const_soln_arrays[box_no]);
        });

    amrex::Gpu::streamSynchronize();
}

void ScalarFieldLevel::specificUpdateODE(amrex::MultiFab &a_soln)
{
    BL_PROFILE("ScalarFieldLevel::specificUpdateODE()");

    const auto &soln_arrays = a_soln.arrays();
    const TraceARemoval trace_A_removal;

    amrex::ParallelFor(a_soln, amrex::IntVect(0),
                       [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
                       { trace_A_removal(ix, iy, iz, soln_arrays[box_no]); });
    amrex::Gpu::streamSynchronize();
}

void ScalarFieldLevel::tag_cells(amrex::TagBoxArray &a_tag_box_array,
                                 const amrex::Real /*a_regrid_threshold*/)
{
    BL_PROFILE("ScalarFieldLevel::tag_cells()");

    const auto &tag_arrays = a_tag_box_array.arrays();
    const FixedGridsTagger tagger(Geom().CellSize(0), Level(),
                                  Geom().ProbLength(0),
                                  simParams().initial_params.center);

    amrex::ParallelFor(a_tag_box_array,
                       [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
                       { tagger(ix, iy, iz, tag_arrays[box_no]); });
    amrex::Gpu::streamSynchronize();
}
