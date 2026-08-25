/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#include "ScalarFieldLevel.hpp"

#include "CCZ4RHSWithMatter.hpp"
#include "ConstraintsWithMatter.hpp"
#include "EMTensor.hpp"
#include "FixedGridsTagger.hpp"
#include "GRParmParse.hpp"
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
using ScalarFieldConstraints =
    ConstraintsWithMatter<ScalarFieldLevel::ScalarFieldWithPotential>;

ScalarFieldAMR *ScalarFieldLevel::get_scalar_field_amr_ptr()
{
    return dynamic_cast<ScalarFieldAMR *>(get_gramr_ptr());
}

void ScalarFieldLevel::variableSetUp()
{
    BL_PROFILE("ScalarFieldLevel::variableSetUp()");
    stateVariableSetUp();
    ScalarFieldConstraints::set_up(state_index);
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

    GRParmParse line_pp("line_extraction");
    bool line_extraction_enabled{};
    line_pp.get("enabled", line_extraction_enabled);

    if (Level() == 0 && line_extraction_enabled)
    {
        const amrex::Real time         = get_state_data(state_index).curTime();
        const amrex::Real dt           = get_gramr_ptr()->dtLevel(0);
        const amrex::Real restart_time = get_gramr_ptr()->get_restart_time();
        const bool first_step          = (time <= dt);

        GRParmParse geometry_pp("geometry");
        std::array<amrex::ParticleReal, AMREX_SPACEDIM> start_coords{};
        geometry_pp.get("center", start_coords);
        auto end_coords = start_coords;
        amrex::ParticleReal max_radius{};
        line_pp.get("max_radius", max_radius);
        const amrex::ParticleReal coordinate_offset =
            max_radius /
            std::sqrt(static_cast<amrex::ParticleReal>(AMREX_SPACEDIM));
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
        {
            end_coords[dir] += coordinate_offset;
        }

        int num_points{};
        line_pp.get("num_points", num_points);
        std::string output_subpath{};
        line_pp.get("output_subpath", output_subpath);
        GRParmParse grteclyn_pp("grteclyn");
        std::string output_path{};
        grteclyn_pp.get("output_path", output_path);
        const std::string data_path = output_path + "/" + output_subpath + "/";

        const LineExtraction<1> phi_extraction(c_phi, num_points, start_coords,
                                               end_coords, dt, time,
                                               restart_time, first_step);
        phi_extraction.execute_query(
            &get_scalar_field_amr_ptr()->phi_interpolator,
            data_path + "phi_profile_");

        const LineExtraction<1> rho_extraction(0, num_points, start_coords,
                                               end_coords, dt, time,
                                               restart_time, first_step);
        const LineExtraction<1>::derived_vars_t rho_vars{
            ScalarFieldEnergyDensity::name, {"rho"}, {BCParity::even}};
        rho_extraction.execute_query(
            &get_scalar_field_amr_ptr()->rho_interpolator,
            data_path + "rho_profile_", rho_vars);
    }
}

void ScalarFieldLevel::initData()
{
    BL_PROFILE("ScalarFieldLevel::initData()");

    if (get_gramr_ptr()->Verbose() > 0)
    {
        amrex::Print() << "ScalarFieldLevel::initData " << Level() << "\n";
    }

    amrex::MultiFab &state_new = get_new_data(state_index);
    const auto &state_arrays   = state_new.arrays();

    OscillatonInitialData::params_t initial_params;
    GRParmParse geometry_pp("geometry");
    geometry_pp.get("center", initial_params.center);
    const OscillatonInitialData initial_data(initial_params,
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
                                       const amrex::Real /*a_time*/)
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

    GRParmParse evolution_pp("evolution");
    int max_spatial_derivative_order{};
    evolution_pp.get("spatial_derivative_order", max_spatial_derivative_order);
    if (max_spatial_derivative_order != 4)
    {
        amrex::Abort("ScalarField currently supports fourth-order spatial "
                     "derivatives only");
    }

    GRParmParse scalar_field_pp("scalar_field");
    Potential::params_t potential_params;
    scalar_field_pp.get("scalar_mass", potential_params.scalar_mass);
    amrex::Real G_Newton{};
    scalar_field_pp.get("G_Newton", G_Newton);
    const Potential potential(potential_params);
    const ScalarFieldWithPotential scalar_field(potential, G_Newton);
    const CCZ4RHSWithMatter<ScalarFieldWithPotential,
                            MovingPunctureGaugeWithMatter,
                            FourthOrderDerivatives>
        ccz4_rhs(scalar_field, Geom().CellSize(0));

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
    GRParmParse geometry_pp("geometry");
    std::array<amrex::Real, AMREX_SPACEDIM> center{};
    geometry_pp.get("center", center);
    const FixedGridsTagger tagger(Geom().CellSize(0), Level(),
                                  Geom().ProbLength(0) * 2.0, center);

    amrex::ParallelFor(a_tag_box_array,
                       [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
                       { tagger(ix, iy, iz, tag_arrays[box_no]); });
    amrex::Gpu::streamSynchronize();
}
