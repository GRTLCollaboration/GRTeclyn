/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#include "BinaryBHLevel.hpp"
#include "BinaryBH.hpp"
#include "CCZ4RHS.hpp"
#include "ChiExtractionTaggingCriterion.hpp"
#include "FilesystemTools.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "PunctureTracker.hpp"
// xxxxx #include "SixthOrderDerivatives.hpp"
#include "TraceARemoval.hpp"
#include "TwoPuncturesInitialData.hpp"
#include "Weyl4.hpp"
#include "WeylExtraction.hpp"

BHAMR *BinaryBHLevel::get_bhamr_ptr()
{
    return dynamic_cast<BHAMR *>(get_gramr_ptr());
}

// Things to do during the advance step after RK4 steps
void BinaryBHLevel::specificAdvance()
{
    amrex::MultiFab &S_new = get_new_data(State_Type);
    const auto &arrs       = S_new.arrays();

    // Enforce the trace free A_ij condition and positive chi and alpha
    amrex::ParallelFor(S_new,
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       {
                           amrex::CellData<amrex::Real> cell =
                               arrs[box_no].cellData(i, j, k);
                           TraceARemoval()(cell);
                           PositiveChiAndAlpha()(cell);
                       });

    // Check for nan's
    if (simParams().nan_check)
    {
        if (S_new.contains_nan(0, S_new.nComp(), amrex::IntVect(0), true))
        {
            amrex::Abort("NaN in specificAdvance");
        }
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
    BinaryBH binary(simParams().bh1_params, simParams().bh2_params,
                    Geom().CellSize(0));

    static_assert(std::is_trivially_copyable_v<BinaryBH>,
                  "BinaryBH needs to be device copyable");

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

    // Enforce positive chi and alpha and trace free A
    amrex::ParallelFor(a_soln, a_soln.nGrowVect(),
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       {
                           amrex::CellData<amrex::Real> cell =
                               soln_arrs[box_no].cellData(i, j, k);
                           TraceARemoval()(cell);
                           PositiveChiAndAlpha()(cell);
                       });

    // Calculate CCZ4 right hand side
    if (simParams().max_spatial_derivative_order == 4)
    {
        CCZ4RHS<MovingPunctureGauge, FourthOrderDerivatives> ccz4rhs(
            simParams().ccz4_params, Geom().CellSize(0), simParams().sigma,
            simParams().formulation);
        amrex::ParallelFor(
            a_rhs,
            [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) {
                ccz4rhs.compute(i, j, k, rhs_arrs[box_no], soln_c_arrs[box_no]);
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
}

void BinaryBHLevel::errorEst(amrex::TagBoxArray &tag_box_array,
                             int /*clearval*/, int /*tagval*/,
                             amrex::Real /*time*/, int /*n_error_buf*/,
                             int /*ngrow*/)
{
    BL_PROFILE("BinaryBHLevel::errorEst()");
    amrex::MultiFab &state_new = get_new_data(State_Type);
    const auto cur_time        = get_state_data(State_Type).curTime();

    const int nghost =
        state_new.nGrow(); // Need ghost cells to compute gradient
    const int ncomp = 1;
    // We only use chi in the tagging criterion so only fill the ghosts for chi
    FillPatch(*this, state_new, nghost, cur_time, State_Type, c_chi, ncomp);

    const auto &simpar = simParams();

    const auto &tags           = tag_box_array.arrays();
    const auto &state_new_arrs = state_new.const_arrays();
    const auto tagval          = amrex::TagBox::SET;

    // TODO: Change to puncture tagging
    ChiExtractionTaggingCriterion tagger(Geom().CellSize(0), Level(),
                                         simpar.extraction_params,
                                         simpar.activate_extraction);
    amrex::Real threshold = simpar.regrid_thresholds[Level()];
    amrex::ParallelFor(state_new, amrex::IntVect(0),
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       {
                           amrex::Real criterion =
                               tagger(i, j, k, state_new_arrs[box_no]);
                           if (criterion >= threshold)
                           {
                               tags[box_no](i, j, k) = tagval;
                           }
                       });
    amrex::Gpu::streamSynchronize();
}

void BinaryBHLevel::restart_punctures()
{
    if (simParams().track_punctures)
    {
        BHAMR *bh_amr_ptr  = get_bhamr_ptr();
        int coarsest_level = 0;
        bh_amr_ptr->m_puncture_tracker.restart(
            bh_amr_ptr->levelSteps(coarsest_level));
    }
}

void BinaryBHLevel::specific_post_init()
{
    BL_PROFILE("BinaryBHLevel::specific_post_init()");

    restart_punctures();
}

void BinaryBHLevel::specific_post_restart()
{
    BL_PROFILE("BinaryBHLevel::specific_post_restart()");

    restart_punctures();
}

void BinaryBHLevel::specificPostCheckpoint(const std::string &a_chk_dir,
                                           std::ostream & /*a_os*/)
{
    if (simParams().track_punctures)
    {
        get_bhamr_ptr()->m_puncture_tracker.checkpoint(a_chk_dir);
    }
}

void BinaryBHLevel::specificPostTimeStep()
{
    // do puncture tracking on requested level
    if (simParams().track_punctures &&
        Level() == simParams().puncture_tracking_level)
    {
        BL_PROFILE("PunctureTracking");
        // only do the write out for every coarsest level timestep
        // int coarsest_level = 0;
        // bool write_punctures = at_level_timestep_multiple(coarsest_level);
        BHAMR *bh_amr            = get_bhamr_ptr();
        bool write_punctures     = true;
        amrex::Real cur_time     = get_state_data(State_Type).curTime();
        amrex::Real restart_time = bh_amr->get_restart_time();
        amrex::Real dt           = bh_amr->dtLevel(Level());
        bh_amr->m_puncture_tracker.execute_tracking(cur_time, restart_time, dt,
                                                    write_punctures);
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
                    VariableType::diagnostic, Interval(c_Weyl4_Re, c_Weyl4_Im),
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
            AMRReductions<VariableType::diagnostic> amr_reductions(m_gr_amr);
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

#ifdef AMREX_USE_HDF5
// Things to do before a plot level - need to calculate the Weyl scalars
void BinaryBHLevel::prePlotLevel()
{
    fillAllGhosts();
    if (m_p.activate_extraction == 1)
    {
        BoxLoops::loop(
            make_compute_pack(
                Weyl4(m_p.extraction_params.center, m_dx, m_p.formulation),
                Constraints(m_dx, c_Ham, Interval(c_Mom1, c_Mom3))),
            m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
    }
}
#endif /* AMREX_USE_HDF5 */
