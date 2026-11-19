/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#include "KerrBHLevel.hpp"

#include "AMRReductions.hpp"
#include "AlgebraicConstraintsEnforcer.hpp"
#include "CCZ4RHS.hpp"
#include "Constraints.hpp"
#include "FilesystemTools.hpp"
#include "FixedGridsTagger.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRParmParse.hpp"
#include "GammaCalculator.hpp"
#include "KerrBHInitialData.hpp"
#include "PositiveChiAndLapse.hpp"
#include "SixthOrderDerivatives.hpp"
#include "SmallDataIO.hpp"

BHAmr<KerrBHLevel::num_punctures> *KerrBHLevel::get_bh_amr_ptr()
{
    return dynamic_cast<BHAmr<num_punctures> *>(get_gr_amr_ptr());
}

void KerrBHLevel::variableSetUp()
{
    BL_PROFILE("KerrBHLevel::variableSetUp()");

    // Set up the state variables
    state_variable_set_up();

    Constraints::set_up(state_index);
}

// Things to do during the advance step after RK4 steps
void KerrBHLevel::specific_advance()
{
    amrex::MultiFab &state_new = get_new_data(state_index);
    const auto &state_arrays   = state_new.arrays();

    // The classes to be used
    AlgebraicConstraintsEnforcer algebraic_constraints_enforcer;
    PositiveChiAndLapse positive_chi_lapse;

    // Enforce det(h)=1, the trace free A_ij condition and positive chi and
    // lapse
    amrex::ParallelFor(state_new,
                       [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
                       {
                           algebraic_constraints_enforcer(ix, iy, iz,
                                                          state_arrays[box_no]);
                           positive_chi_lapse(ix, iy, iz, state_arrays[box_no]);
                       });
}

// Initial Data setup for the Kerr Black Hole
void KerrBHLevel::initData()
{
    BL_PROFILE("KerrBHLevel::initData");
    if (get_gr_amr_ptr()->Verbose() > 0)
    {
        amrex::Print() << "KerrBHLevel::initData " << Level() << "\n";
    }

    // Set up the compute class for the KerrBH initial data
    amrex::Real dx = Geom().CellSize(0);
    KerrBHInitialData kerr_initial_data(dx);
    static_assert(std::is_trivially_copyable_v<KerrBHInitialData>,
                  "KerrBHInitialData needs to be device copyable");

    // First set everything to zero (to avoid undefinded values in constraints)
    // then calculate initial data
    amrex::MultiFab &state_new = get_new_data(state_index);
    const auto &state_arrays   = state_new.arrays();
    amrex::ParallelFor(state_new, state_new.nGrowVect(),
                       [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
                       {
                           amrex::CellData<amrex::Real> cell =
                               state_arrays[box_no].cellData(ix, iy, iz);
                           for (int n = 0; n < cell.nComp(); ++n)
                           {
                               cell[n] = 0.;
                           }
                           kerr_initial_data(ix, iy, iz, state_arrays[box_no]);
                       });

    // Gamma^i needs derivatives of the metric, so it is calculated after the
    // initial data has been set everywhere (including the ghost cells)
    if (m_evolution_spatial_derivative_order == 4)
    {
        const GammaCalculator<FourthOrderDerivatives> gamma_calculator(
            Geom().CellSize(0));
        amrex::ParallelFor(
            state_new, [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
            { gamma_calculator(ix, iy, iz, state_arrays[box_no]); });
    }
    else if (m_evolution_spatial_derivative_order == 6)
    {
        const GammaCalculator<SixthOrderDerivatives> gamma_calculator(
            Geom().CellSize(0));
        amrex::ParallelFor(
            state_new, [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
            { gamma_calculator(ix, iy, iz, state_arrays[box_no]); });
    }

    amrex::Gpu::streamSynchronize();
}

// Calculate RHS during RK4 substeps
// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
void KerrBHLevel::specific_eval_rhs(amrex::MultiFab &a_soln,
                                    amrex::MultiFab &a_rhs,
                                    const amrex::Real /*a_time*/)
{
    BL_PROFILE("KerrBHLevel::specific_eval_rhs()");
    const auto &soln_arrays       = a_soln.arrays();
    const auto &const_soln_arrays = a_soln.const_arrays();
    const auto &rhs_arrays        = a_rhs.arrays();
    const auto soln_ghosts        = a_soln.nGrowVect();

    // The classes to be used
    AlgebraicConstraintsEnforcer algebraic_constraints_enforcer;
    PositiveChiAndLapse positive_chi_lapse;

    // Enforce positive chi and lapse, det(h)=1 and trace free A
    amrex::ParallelFor(a_soln, soln_ghosts,
                       [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
                       {
                           algebraic_constraints_enforcer(ix, iy, iz,
                                                          soln_arrays[box_no]);
                           positive_chi_lapse(ix, iy, iz, soln_arrays[box_no]);
                       });

    // Calculate CCZ4 right hand side using dynamic derivative order
    if (m_evolution_spatial_derivative_order == 4)
    {
        CCZ4RHS<FourthOrderDerivatives> ccz4rhs(Geom().CellSize(0));
        MovingPunctureGauge<FourthOrderDerivatives> moving_puncture_gauge(
            Geom().CellSize(0));

        // NB: These are split up to avoid having to pre-compute all the
        //  first and second derivatives in memory on the GPU at once.

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
    }
    else if (m_evolution_spatial_derivative_order == 6)
    {
        CCZ4RHS<SixthOrderDerivatives> ccz4rhs(Geom().CellSize(0));
        MovingPunctureGauge<SixthOrderDerivatives> moving_puncture_gauge(
            Geom().CellSize(0));

        // NB: These are split up to avoid having to pre-compute all the
        //  first and second derivatives in memory on the GPU at once.

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
    }

    amrex::Gpu::streamSynchronize();
}

// enforce algebraic constraints during RK4 substeps
void KerrBHLevel::specific_update_ode(amrex::MultiFab &a_soln)
{

    AlgebraicConstraintsEnforcer algebraic_constraints_enforcer;
    const auto soln_ghosts = amrex::IntVect(0); // zero ghost cells

    // Enforce the det(h)=1 and trace free A_ij conditions
    const auto &soln_arrays = a_soln.arrays();
    amrex::ParallelFor(
        a_soln, soln_ghosts,
        [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
        { algebraic_constraints_enforcer(ix, iy, iz, soln_arrays[box_no]); });

    amrex::Gpu::streamSynchronize();
}

void KerrBHLevel::tag_cells(amrex::TagBoxArray &a_tag_box_array,
                            amrex::Real /*a_regrid_threshold*/)
{
    BL_PROFILE("KerrBHLevel::tag_cells()");

    const auto &tag_arrays = a_tag_box_array.arrays();

    const FixedGridsTagger tagger(Geom().CellSize(0), Level());

    amrex::ParallelFor(a_tag_box_array,
                       [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
                       { tagger(ix, iy, iz, tag_arrays[box_no]); });
    amrex::Gpu::streamSynchronize();
}

void KerrBHLevel::specific_post_timestep()
{
    BL_PROFILE("KerrBHLevel::specific_post_timestep");

    if (Level() == 0)
    {
        GRParmParse pp;
        bool store_constraint_norms = false;
        pp.query("diagnostic.store_constraint_norms", store_constraint_norms);
        if (store_constraint_norms)
        {
            write_constraint_norms();
        }
    }
}

void KerrBHLevel::write_constraint_norms()
{
    BL_PROFILE("KerrBHLevel::write_constraint_norms()");

    const amrex::Real time         = get_state_data(state_index).curTime();
    const amrex::Real dt           = get_gr_amr_ptr()->dtLevel(0);
    const amrex::Real restart_time = get_gr_amr_ptr()->get_restart_time();
    const bool first_step          = (time <= dt);

    // Only one derived quantity is reduced, so it can be derived on all levels
    // directly (to reduce several quantities in one pass over the grid, derive
    // them together with an AMRReductions::DerivedHierarchy instead).
    // Component 0 is Ham and components 1 to AMREX_SPACEDIM are Mom1-3 (see
    // Constraints::var_names); the latter are combined into a single Mom norm.
    auto &amr              = *get_gr_amr_ptr();
    const auto constraints = amr.derive(Constraints::name, time, 0);

    const AMRReductions::AMRReducer reducer(amr);
    const auto [ham_l2_norm, mom_l2_norm] = reducer.compute(
        amrex::GetVecOfConstPtrs(constraints), AMRReductions::L2Norm{.comp = 0},
        AMRReductions::L2Norm{
            .comp      = 1,
            .num_comps = AMREX_SPACEDIM,
        });

    GRParmParse pp;
    std::string output_path;
    pp.get("grteclyn.output_path", output_path);
    const std::string diagnostics_path = output_path + "/diagnostics";
    FilesystemTools::ensure_directory_exists(diagnostics_path);

    SmallDataIO constraints_file(diagnostics_path + "/constraint_norms", dt,
                                 time, restart_time, SmallDataIO::APPEND,
                                 first_step);
    if (first_step)
    {
        constraints_file.write_header_line({"L2_Ham", "L2_Mom"});
    }
    constraints_file.write_time_data_line(
        std::vector<amrex::Real>{ham_l2_norm, mom_l2_norm});
}
