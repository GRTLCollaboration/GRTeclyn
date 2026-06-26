/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#include "BinaryBHLevel.hpp"
#include "BinaryBHInitialData.hpp"
#include "CCZ4RHS.hpp"
#include "ChiTagger.hpp"
#include "Constraints.hpp"
#include "ExtractionTagger.hpp"
#include "GaugeFixer.hpp"
#include "PositiveChiAndLapse.hpp"
#include "PunctureTagger.hpp"
#include "PunctureTracker.hpp"
#include "SixthOrderDerivatives.hpp"
#include "TraceARemoval.hpp"
#include "TwoPuncturesInitialData.hpp"
#include "Weyl4.hpp"
#include "WeylExtraction.hpp"

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

    Constraints::set_up(state_index);

    Weyl4::set_up(state_index);
}

// Things to do during the advance step after RK4 steps
void BinaryBHLevel::specificAdvance()
{
    amrex::MultiFab &state_new = get_new_data(state_index);
    const auto &state_arrays   = state_new.arrays();

    // The classes to be used
    TraceARemoval trace_A_removal;
    PositiveChiAndLapse positive_chi_lapse;

    // Enforce the trace free A_ij condition and positive chi and lapse
    amrex::ParallelFor(state_new,
                       [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
                       {
                           trace_A_removal(ix, iy, iz, state_arrays[box_no]);
                           positive_chi_lapse(ix, iy, iz, state_arrays[box_no]);
                       });
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
    TwoPuncturesInitialData two_punctures_initial_data(Geom().CellSize(0),
                                                       simParams().center);

    two_punctures_initial_data.solve(); // only solves first time

    // hack to modify const parameters before we refactor parameters
    // auto sim_params = simParams();
    // two_punctures_initial_data.set_bh_params(sim_params.bh1_params,
    //                                          sim_params.bh2_params);
    // get_gramr_ptr()->set_simulation_parameters(sim_params);

    amrex::MultiFab &state_new = get_new_data(state_index);
#ifdef AMREX_USE_GPU
    amrex::MFInfo mf_info;
    mf_info.SetArena(amrex::The_Cpu_Arena());
    amrex::MultiFab host_state(state_new.boxArray(),
                               state_new.DistributionMap(), state_new.nComp(),
                               state_new.nGrowVect(), mf_info);
#else
    amrex::MultiFab &host_state = state_new;
#endif

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(state_new, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi)
    {
        const amrex::Box &grown_tile_box = mfi.growntilebox();
        const auto &state_array          = host_state.array(mfi);

        amrex::LoopOnCpu(
            grown_tile_box, [=](int ix, int iy, int iz)
            { two_punctures_initial_data(ix, iy, iz, state_array); });
#ifdef AMREX_USE_GPU
        // Copy to device
        amrex::Gpu::htod_memcpy_async(
            state_new[mfi].dataPtr(), host_state[mfi].dataPtr(),
            host_state[mfi].size() * sizeof(amrex::Real));
#endif
    }

#else
    // Set up the compute class for the BinaryBH initial data
    double dx = Geom().CellSize(0);
    BinaryBHInitialData binary_initial_data(simParams().bh1_params,
                                            simParams().bh2_params, dx);

    static_assert(std::is_trivially_copyable_v<BinaryBHInitialData>,
                  "BinaryBHInitialData needs to be device copyable");

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
                           binary_initial_data(ix, iy, iz,
                                               state_arrays[box_no]);
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
// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
void BinaryBHLevel::specificEvalRHS(amrex::MultiFab &a_soln,
                                    amrex::MultiFab &a_rhs,
                                    const double /*a_time*/)
{
    BL_PROFILE("BinaryBHLevel::specificEvalRHS()");
    const auto &soln_arrays       = a_soln.arrays();
    const auto &const_soln_arrays = a_soln.const_arrays();
    const auto &rhs_arrays        = a_rhs.arrays();
    const auto soln_ghosts        = a_soln.nGrowVect();

    // The classes to be used
    TraceARemoval trace_A_removal;
    PositiveChiAndLapse positive_chi_lapse;

    // Enforce positive chi and lapse and trace free A
    amrex::ParallelFor(a_soln, soln_ghosts,
                       [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
                       {
                           trace_A_removal(ix, iy, iz, soln_arrays[box_no]);
                           positive_chi_lapse(ix, iy, iz, soln_arrays[box_no]);
                       });

    // Calculate CCZ4 right hand side
    if (simParams().max_spatial_derivative_order == 4)
    {
        CCZ4RHS<MovingPunctureGauge<FourthOrderDerivatives>, FourthOrderDerivatives> ccz4rhs(
            simParams().ccz4_params, Geom().CellSize(0), simParams().sigma,
            simParams().formulation);

        amrex::ParallelFor(
            a_rhs,
            [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
            {
                ccz4rhs(ix, iy, iz, rhs_arrays[box_no],
                        const_soln_arrays[box_no]);
            });
    }
    else if (simParams().max_spatial_derivative_order == 6)
    {
        CCZ4RHS<MovingPunctureGauge<SixthOrderDerivatives>, SixthOrderDerivatives> ccz4rhs(
            simParams().ccz4_params_d, Geom().CellSize(0), simParams().sigma,
            simParams().formulation);

        amrex::ParallelFor(
            a_rhs,
            [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
            {
                ccz4rhs(ix, iy, iz, rhs_arrays[box_no],
                        const_soln_arrays[box_no]);
            });
    }

    GaugeFixer gaugefix(Geom().CellSize(0), simParams().center);

    amrex::ParallelFor(
            a_rhs,
            [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
            {
                gaugefix(ix, iy, iz, rhs_arrays[box_no]);
            });

    amrex::Gpu::streamSynchronize();


}

// enforce trace removal during RK4 substeps
void BinaryBHLevel::specificUpdateODE(amrex::MultiFab &a_soln)
{

    TraceARemoval trace_A_removal;
    const auto soln_ghosts = amrex::IntVect(0); // zero ghost cells

    // Enforce the trace free A_ij condition
    const auto &soln_arrays = a_soln.arrays();
    amrex::ParallelFor(a_soln, soln_ghosts,
                       [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
                       { trace_A_removal(ix, iy, iz, soln_arrays[box_no]); });

    amrex::Gpu::streamSynchronize();
}

void BinaryBHLevel::pre_tag_cells()
{
    amrex::MultiFab &state_new = get_new_data(state_index);
    const auto current_time    = get_state_data(state_index).curTime();

    // Just fill 2 ghosts for chi to calculate second derivatives
    const int nghost = 2;
    const int ncomp  = 1;
    FillPatch(*this, state_new, nghost, current_time, state_index, c_chi,
              ncomp);
}

void BinaryBHLevel::tag_cells(amrex::TagBoxArray &a_tag_box_array,
                              amrex::Real a_regrid_threshold)
{
    BL_PROFILE("BinaryBHLevel::tag_cells()");
    amrex::MultiFab &state_new = get_new_data(state_index);

    const auto &tag_arrays         = a_tag_box_array.arrays();
    const auto &state_const_arrays = state_new.const_arrays();

    ChiTagger chi_tagger(Geom().CellSize(0), a_regrid_threshold);

    ExtractionTagger extraction_tagger(Geom().CellSize(0), Level(),
                                       simParams().extraction_params,
                                       simParams().activate_extraction);

    const bool puncture_tracking_enabled =
        simParams().puncture_tracking_enabled;
    constexpr auto num_puncture_coords =
        static_cast<std::size_t>(AMREX_SPACEDIM * num_punctures);
    std::array<amrex::Real, num_puncture_coords> puncture_coords{};

    if (puncture_tracking_enabled)
    {
        puncture_coords = get_puncture_tracker().get_puncture_coords();
    }

    PunctureTagger<num_punctures> puncture_tagger(
        Geom().CellSize(0), Level(), get_gramr_ptr()->maxLevel(),
        puncture_coords,
        {simParams().bh1_params.mass, simParams().bh2_params.mass});

    amrex::ParallelFor(state_new, amrex::IntVect(0),
                       [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
                       {
                           chi_tagger(ix, iy, iz, tag_arrays[box_no],
                                      state_const_arrays[box_no]);

                           extraction_tagger(ix, iy, iz, tag_arrays[box_no]);

                           if (puncture_tracking_enabled)
                           {
                               puncture_tagger(ix, iy, iz, tag_arrays[box_no]);
                           }
                       });

    amrex::Gpu::streamSynchronize();
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
        amrex::Real current_time = get_state_data(state_index).curTime();
        amrex::Real dt           = get_gramr_ptr()->dtLevel(Level());
        get_puncture_tracker().track(current_time, dt, write_punctures);
    }

    // Weyl extraction
    if (simParams().activate_extraction)
    {
        int min_level = simParams().extraction_params.min_extraction_level();
        bool calculate_weyl = at_level_timestep_multiple(min_level);

        if (calculate_weyl && Level() == min_level)
        {
            amrex::Real m_time       = get_state_data(state_index).curTime();
            amrex::Real m_dt         = get_gramr_ptr()->dtLevel(Level());
            amrex::Real restart_time = get_gramr_ptr()->get_restart_time();
            bool first_step          = (m_time <= m_dt);

            WeylExtraction my_extraction(simParams().extraction_params, m_dt,
                                         m_time, first_step, restart_time);
            my_extraction.execute_query(&get_bhamr_ptr()->m_weyl_interpolator,
                                        "Weyl4");
        }
    }

#if 0
//xxxxx specificPostTimeStep
    BL_PROFILE("BinaryBHLevel::specificPostTimeStep");

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
