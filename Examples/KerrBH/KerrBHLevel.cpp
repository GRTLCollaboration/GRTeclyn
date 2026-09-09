/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#include "KerrBHLevel.hpp"

#include "AlgebraicConstraintsEnforcer.hpp"
#include "KerrBHInitialData.hpp"
#include "CCZ4RHS.hpp"
#include "ChiTagger.hpp"
#include "Constraints.hpp"
#include "ExtractionTagger.hpp"
#include "FourthOrderDerivatives.hpp"
#include "PositiveChiAndLapse.hpp"
#include "SixthOrderDerivatives.hpp"
#include "Weyl4.hpp"
#include "WeylExtraction.hpp"

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

    Weyl4::set_up(state_index);
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
                           kerr_initial_data(ix, iy, iz,
                                               state_arrays[box_no]);
                       });

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

void KerrBHLevel::pre_tag_cells()
{
    amrex::MultiFab &state_new = get_new_data(state_index);
    const auto current_time    = get_state_data(state_index).curTime();

    // Fill ghosts for chi to calculate second derivatives
    // 4th-order d2 requires 2 ghost cells
    const int num_ghosts = 2;
    const int num_comps  = 1;

    FillPatch(*this, state_new, num_ghosts, current_time, state_index, c_chi,
              num_comps);
}

void KerrBHLevel::tag_cells(amrex::TagBoxArray &a_tag_box_array,
                              amrex::Real a_regrid_threshold)
{
    BL_PROFILE("KerrBHLevel::tag_cells()");
    amrex::MultiFab &state_new = get_new_data(state_index);

    const auto &tag_arrays         = a_tag_box_array.arrays();
    const auto &state_const_arrays = state_new.const_arrays();

    ChiTagger chi_tagger(Geom().CellSize(0), a_regrid_threshold);

    spherical_extraction_params_t extraction_params("weyl_extraction");
    extraction_params.fill_params();
    ExtractionTagger extraction_tagger(Geom().CellSize(0), Level(),
                                       extraction_params);

    amrex::ParallelFor(state_new, amrex::IntVect(0),
                       [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
                       {
                           chi_tagger(ix, iy, iz, tag_arrays[box_no],
                                      state_const_arrays[box_no]);

                           extraction_tagger(ix, iy, iz, tag_arrays[box_no]);
                       });

    amrex::Gpu::streamSynchronize();
}

void KerrBHLevel::specific_post_init()
{
    BL_PROFILE("KerrBHLevel::specific_post_init()");

}

void KerrBHLevel::specific_post_restart()
{
    BL_PROFILE("KerrBHLevel::specific_post_restart()");

}

void KerrBHLevel::specific_post_plotfile(const std::string &a_dir,
                                           std::ostream &a_os)
{
    BL_PROFILE("KerrBHLevel::specific_post_plotfile()");
}

void KerrBHLevel::specific_post_checkpoint(const std::string &a_chk_dir,
                                             std::ostream & /*a_os*/)
{
    BL_PROFILE("KerrBHLevel::specific_post_checkpoint()");
}

void KerrBHLevel::specific_post_timestep()
{
    BL_PROFILE("KerrBHLevel::specific_post_timestep");

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
            my_extraction.execute_query(&get_bh_amr_ptr()->m_weyl_interpolator);
        }
    }
}
