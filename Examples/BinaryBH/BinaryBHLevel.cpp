/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#include "BinaryBHLevel.hpp"
#include "BinaryBHInitialData.hpp"
#include "CCZ4RHS.hpp"
#include "ChiExtractionTagger.hpp"
#include "Constraints.hpp"
#include "PositiveChiAndLapse.hpp"
#include "PunctureTagger.hpp"
#include "PunctureTracker.hpp"
// xxxxx #include "SixthOrderDerivatives.hpp"
#include "TraceARemoval.hpp"
#include "TwoPuncturesInitialData.hpp"
#include "Weyl4.hpp"

BHAMR<BinaryBHLevel::num_punctures> *BinaryBHLevel::get_bhamr_ptr()
{
    return dynamic_cast<BHAMR<num_punctures> *>(get_gramr_ptr());
}

PunctureTracker<BinaryBHLevel::num_punctures> &
BinaryBHLevel::get_puncture_tracker()
{
    return get_bhamr_ptr()->get_puncture_tracker();
}

void BinaryBHLevel::variableSetUp()
{
    BL_PROFILE("BinaryBHLevel::variableSetUp()");

    // Set up the state variables
    stateVariableSetUp();

    Constraints::set_up(State_Type);

    Weyl4::set_up(State_Type);
}

// Things to do during the advance step after RK4 steps
void BinaryBHLevel::specificAdvance()
{
    amrex::MultiFab &S_new = get_new_data(State_Type);
    const auto &arrs       = S_new.arrays();

    // Enforce the trace free A_ij condition and positive chi and lapse
    amrex::ParallelFor(S_new,
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       {
                           amrex::CellData<amrex::Real> cell =
                               arrs[box_no].cellData(i, j, k);
                           TraceARemoval()(cell);
                           PositiveChiAndLapse()(cell);
                       });

    // Check for nan's
    if (simParams().nan_check)
    {
        if (S_new.contains_nan(0, S_new.nComp(), amrex::IntVect(0), true))
        {
            amrex::Abort("NaN in specificAdvance");
        }
    }
    else
    {
        // stream sync already present in nan_check so only need this here
        amrex::Gpu::streamSynchronize();
    }
}

// This initial data uses an approximation for the metric which
// is valid for small boosts
void BinaryBHLevel::initData()
{
    BL_PROFILE("BinaryBHLevel::initialData");
    if (m_verbosity > 0)
    {
        amrex::Print() << "BinaryBHLevel::initialData " << Level() << "\n";
    }
#ifdef USE_TWOPUNCTURES
    // xxxxx USE_TWOPUNCTURES todo
    TwoPuncturesInitialData two_punctures_initial_data(
        m_dx, m_p.center, m_tp_amr.m_two_punctures);
    // Can't use simd with this initial data
    BoxLoops::loop(two_punctures_initial_data, m_state_new, m_state_new,
                   INCLUDE_GHOST_CELLS, disable_simd());
#else
    // Set up the compute class for the BinaryBH initial data
    BinaryBHInitialData binary(simParams().bh1_params, simParams().bh2_params,
                               Geom().CellSize(0));

    static_assert(std::is_trivially_copyable_v<BinaryBHInitialData>,
                  "BinaryBHInitialData needs to be device copyable");

    // First set everything to zero (to avoid undefinded values in constraints)
    // then calculate initial data
    amrex::MultiFab &state = get_new_data(State_Type);
    const auto &arrs       = state.arrays();
    amrex::ParallelFor(state, state.nGrowVect(),
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       {
                           amrex::CellData<amrex::Real> cell =
                               arrs[box_no].cellData(i, j, k);
                           for (int n = 0; n < cell.nComp(); ++n)
                           {
                               cell[n] = 0.;
                           }
                           binary.init_data(i, j, k, cell);
                       });
#endif
    amrex::Gpu::streamSynchronize();

    if (simParams().puncture_tracking_enabled && Level() == 0)
    {
        // need to set the puncture coordinates as we use it for the puncture
        // tagging
        get_puncture_tracker().set_puncture_coords(
            {simParams().bh1_params.center[0], simParams().bh1_params.center[1],
             simParams().bh1_params.center[2], simParams().bh2_params.center[0],
             simParams().bh2_params.center[1],
             simParams().bh2_params.center[2]});
        // can't call start_from_initial_punctures() because we need the full
        // AMR grid first
    }
}

// Calculate RHS during RK4 substeps
void BinaryBHLevel::specificEvalRHS(amrex::MultiFab &a_soln,
                                    amrex::MultiFab &a_rhs,
                                    const double /*a_time*/)
{
    BL_PROFILE("BinaryBHLevel::specificEvalRHS()");
    const auto &soln_arrs   = a_soln.arrays();
    const auto &soln_c_arrs = a_soln.const_arrays();
    const auto &rhs_arrs    = a_rhs.arrays();

    // Enforce positive chi and lapse and trace free A
    amrex::ParallelFor(a_soln, a_soln.nGrowVect(),
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       {
                           amrex::CellData<amrex::Real> cell =
                               soln_arrs[box_no].cellData(i, j, k);
                           TraceARemoval()(cell);
                           PositiveChiAndLapse()(cell);
                       });

    // Calculate CCZ4 right hand side
    if (simParams().max_spatial_derivative_order == 4)
    {
        CCZ4RHS<MovingPunctureGauge, FourthOrderDerivatives> ccz4rhs(
            simParams().ccz4_params, Geom().CellSize(0), simParams().sigma,
            simParams().formulation);
        amrex::ParallelFor(a_rhs,
                           [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                           {
                               ccz4rhs.compute(i, j, k, rhs_arrs[box_no],
                                               soln_c_arrs[box_no]);
                           });
    }
    else if (simParams().max_spatial_derivative_order == 6)
    {
        amrex::Abort("xxxxx max_spatial_derivative_order == 6 todo");
#if 0
        CCZ4RHS<MovingPunctureGauge, SixthOrderDerivatives>
            ccz4rhs(simParams().ccz4_params, Geom().CellSize(0), simParams().sigma,
                    simParams().formulation);
        amrex::ParallelFor(a_rhs,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k)
        {
            amrex::CellData<amrex::Real const> state = soln_c_arrs[box_no].cellData(i,j,k);
            amrex::CellData<amrex::Real> rhs = rhs_arrs[box_no].cellData(i,j,k);
            ccz4rhs.compute(rhs, state);
        });
#endif
    }

    amrex::Gpu::streamSynchronize();
}

// enforce trace removal during RK4 substeps
void BinaryBHLevel::specificUpdateODE(amrex::MultiFab &a_soln)
{
    // Enforce the trace free A_ij condition
    const auto &soln_arrs = a_soln.arrays();
    amrex::ParallelFor(a_soln, amrex::IntVect(0), // zero ghost cells
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       {
                           amrex::CellData<amrex::Real> cell =
                               soln_arrs[box_no].cellData(i, j, k);
                           TraceARemoval()(cell);
                       });

    amrex::Gpu::streamSynchronize();
}

void BinaryBHLevel::pre_tag_cells()
{
    amrex::MultiFab &state_new = get_new_data(State_Type);
    const auto cur_time        = get_state_data(State_Type).curTime();

    // Just fill 2 ghosts for chi to calculate second derivatives
    const int nghost = 2;
    const int ncomp  = 1;
    FillPatch(*this, state_new, nghost, cur_time, State_Type, c_chi, ncomp);
}

void BinaryBHLevel::tag_cells(amrex::TagBoxArray &a_tag_box_array,
                              amrex::Real a_regrid_threshold)
{
    BL_PROFILE("BinaryBHLevel::tag_cells()");
    amrex::MultiFab &state_new = get_new_data(State_Type);

    const auto &tag_arrs       = a_tag_box_array.arrays();
    const auto &state_new_arrs = state_new.const_arrays();

    ChiExtractionTagger chi_extraction_tagger(
        Geom().CellSize(0), Level(), a_regrid_threshold,
        simParams().extraction_params, simParams().activate_extraction);

    const bool puncture_tracking_enabled =
        simParams().puncture_tracking_enabled;
    constexpr auto num_puncture_coords =
        static_cast<std::size_t>(AMREX_SPACEDIM * num_punctures);
    std::array<amrex::Real, num_puncture_coords> puncture_coords{};

    if (puncture_tracking_enabled)
    {
        puncture_coords = get_puncture_tracker().get_puncture_coords();
    }

    // Even though we create this object, it won't be used if puncture tracking
    // is not enabled.
    PunctureTagger<num_punctures> puncture_tagger(
        Geom().CellSize(0), Level(), get_gramr_ptr()->maxLevel(),
        puncture_coords,
        {simParams().bh1_params.mass, simParams().bh2_params.mass});

    amrex::ParallelFor(state_new, amrex::IntVect(0),
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       {
                           chi_extraction_tagger(i, j, k, tag_arrs[box_no],
                                                 state_new_arrs[box_no]);
                           if (puncture_tracking_enabled)
                           {
                               puncture_tagger(i, j, k, tag_arrs[box_no]);
                           }
                       });
    amrex::Gpu::streamSynchronize();
}

// This is called by amrex on every level after regridding. After regrid we
// should redistribute the particles again if BoxArrays or DistributionMappings
// have changed, for example?
void BinaryBHLevel::specific_post_regrid(int a_lbase, int a_new_finest)
{
    // amrex::Print()
    //     << "BinaryBHLevel::specific_post_regrid() on level " << Level()
    //     << "\n";

    // NOTE: this will be called on each level, so redistribute() will
    // happen several times. This is lossy, should be ideally
    // changed so that it happens only once.

    auto *gramr = get_gramr_ptr();
    auto *bh    = get_bhamr_ptr();

    // amrex::Print() << "cumTime = " << (gramr ? gramr->cumTime() : -1.0) <<
    // "\n"; amrex::Print() << "bh ptr  = " << (void *)bh << "\n";

    if (bh && bh->m_weyl_interpolator)
    {
        // amrex::Print() << "Forcing redistribute flag on this rank/level\n";
        bh->m_weyl_interpolator->force_redistribute(true);
        bh->m_weyl_interpolator->Redistribute();
    }
    // else
    // {
    //     amrex::Print() << "Skipping: interpolator is null here\n";
    // }
}

void BinaryBHLevel::specific_post_init()
{
    BL_PROFILE("BinaryBHLevel::specific_post_init()");

    if (simParams().puncture_tracking_enabled)
    {
        get_puncture_tracker().start_from_initial_punctures();
    }
}

void BinaryBHLevel::specific_post_restart()
{
    BL_PROFILE("BinaryBHLevel::specific_post_restart()");

    if (simParams().puncture_tracking_enabled)
    {
        std::string restart_checkpoint{};
        GRParmParse pp("amr");
        pp.get("restart", restart_checkpoint);
        get_puncture_tracker().restart(restart_checkpoint);
    }
}

void BinaryBHLevel::specific_post_plotfile(const std::string &a_dir,
                                           std::ostream &a_os)
{
    if (simParams().puncture_tracking_enabled)
    {
        get_puncture_tracker().write_plotfile(a_dir);
    }
}

void BinaryBHLevel::specific_post_checkpoint(const std::string &a_chk_dir,
                                             std::ostream & /*a_os*/)
{
    if (simParams().puncture_tracking_enabled)
    {
        get_puncture_tracker().checkpoint(a_chk_dir);
    }
}

void BinaryBHLevel::specificPostTimeStep()
{

    // do puncture tracking on requested level
    if (simParams().puncture_tracking_enabled &&
        Level() == simParams().puncture_tracking_level)
    {
        BL_PROFILE("PunctureTracking");

        // only do the write out when we're at at a multiple of the
        // writeout_level
        bool write_punctures = at_level_timestep_multiple(
            simParams().puncture_tracking_writeout_level);
        amrex::Real cur_time = get_state_data(State_Type).curTime();
        amrex::Real dt       = get_gramr_ptr()->dtLevel(Level());
        get_puncture_tracker().track(cur_time, dt, write_punctures);
    }

#if 0
//xxxxx specificPostTimeStep
    BL_PROFILE("BinaryBHLevel::specificPostTimeStep");

    bool first_step =
        (m_time == 0.); // this form is used when 'specificPostTimeStep' was
                        // called during setup at t=0 from Main
    // bool first_step = (m_time == m_dt); // if not called in Main

    if (m_p.activate_extraction == 1)
    {
        int min_level = m_p.extraction_params.min_extraction_level();
        bool calculate_weyl = at_level_timestep_multiple(min_level);
        if (calculate_weyl)
        {
            // Populate the Weyl Scalar values on the grid
            fillAllGhosts();
            BoxLoops::loop(
                Weyl4(m_p.extraction_params.center, m_dx, m_p.formulation),
                m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);

            // Do the extraction on the min extraction level
            if (m_level == min_level)
            {
                BL_PROFILE("WeylExtraction");
                // Now refresh the interpolator and do the interpolation
                // fill ghosts manually to minimise communication
                bool fill_ghosts = false;
                m_gr_amr.m_interpolator->refresh(fill_ghosts);
                m_gr_amr.fill_multilevel_ghosts(
                    VariableType::derived, Interval(c_Weyl4_Re, c_Weyl4_Im),
                    min_level);
                WeylExtraction my_extraction(m_p.extraction_params, m_dt,
                                             m_time, first_step,
                                             m_restart_time);
                my_extraction.execute_query(m_gr_amr.m_interpolator);
            }
        }
    }

    if (m_p.calculate_constraint_norms)
    {
        fillAllGhosts();
        BoxLoops::loop(Constraints(m_dx, c_Ham, Interval(c_Mom1, c_Mom3)),
                       m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
        if (m_level == 0)
        {
            AMRReductions<VariableType::derived> amr_reductions(m_gr_amr);
            double L2_Ham = amr_reductions.norm(c_Ham);
            double L2_Mom = amr_reductions.norm(Interval(c_Mom1, c_Mom3));
            SmallDataIO constraints_file(m_p.data_path + "constraint_norms",
                                         m_dt, m_time, m_restart_time,
                                         SmallDataIO::APPEND, first_step);
            constraints_file.remove_duplicate_time_data();
            if (first_step)
            {
                constraints_file.write_header_line({"L^2_Ham", "L^2_Mom"});
            }
            constraints_file.write_time_data_line({L2_Ham, L2_Mom});
        }
    }

#endif
}
