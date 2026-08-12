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
#include "PositiveChiAndLapse.hpp"
#include "PunctureTagger.hpp"
#include "PunctureTracker.hpp"
// xxxxx #include "SixthOrderDerivatives.hpp"
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
    if (get_gramr_ptr()->Verbose() > 0)
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
    double dx = Geom().CellSize(0);
    BinaryBHInitialData binary_initial_data(dx);
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

    if (get_bhamr_ptr()->puncture_tracking_enabled && Level() == 0)
    {
        // need to set the puncture coordinates as we use it for the puncture
        // tagging
        GRParmParse pp;
        std::array<double, AMREX_SPACEDIM> bh1_center;
        std::array<double, AMREX_SPACEDIM> bh2_center;
        pp.get("bh1.center", bh1_center);
        pp.get("bh2.center", bh2_center);

        get_puncture_tracker().set_puncture_coords(
            {bh1_center[0], bh1_center[1], bh1_center[2], bh2_center[0],
             bh2_center[1], bh2_center[2]});
        // can't call start_from_initial_punctures() because we need the full
        // AMR grid first
    }
}

// Calculate RHS during RK4 substeps
void BinaryBHLevel::specificEvalRHS(amrex::MultiFab &a_soln,
                                    amrex::MultiFab &a_rhs,
                                    const double /*a_time*/)
{
    GRParmParse pp;

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
    int max_spatial_derivative_order;
    pp.get("amr.max_spatial_derivative_order", max_spatial_derivative_order);

    if (max_spatial_derivative_order == 4)
    {
        CCZ4RHS<MovingPunctureGauge, FourthOrderDerivatives> ccz4rhs(
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

        enum formulations : int
        {
            USE_CCZ4 = 0,
            USE_BSSN = 1
        };

        enum covariantZ4 : int
        {
            YES,
            NO
        };

        int use_bssn{0};
        int use_covariantZ4{1};
        pp.query("formulation", use_bssn);
        pp.query("covariantZ4", use_covariantZ4);

        amrex::AnyCTO(
            amrex::TypeList<
                amrex::CompileTimeOptions<USE_CCZ4, USE_BSSN>,
                amrex::CompileTimeOptions<covariantZ4::YES, covariantZ4::NO>>{},
            {use_bssn, use_covariantZ4},
            [&](auto cto_func) { amrex::ParallelFor(a_rhs, cto_func); },
            [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz,
                                 auto formulation, auto control)
            {
                //
                ccz4rhs.compute_A_ij_and_Theta_and_Gamma<formulation, control>(
                    ix, iy, iz, rhs_arrays[box_no], const_soln_arrays[box_no]);
            });

        amrex::ParallelFor(
            a_rhs,
            [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
            {
                ccz4rhs.apply_gauge(ix, iy, iz, rhs_arrays[box_no],
                                    const_soln_arrays[box_no]);

                ccz4rhs.apply_dissipation(ix, iy, iz, rhs_arrays[box_no],
                                          const_soln_arrays[box_no]);
            });
    }
    else if (max_spatial_derivative_order == 6)
    {
        amrex::Abort("xxxxx max_spatial_derivative_order == 6 todo");
#if 0
        CCZ4RHS<MovingPunctureGauge, SixthOrderDerivatives>
            ccz4rhs(Geom().CellSize(0));
        amrex::ParallelFor(a_rhs,
        [=] AMREX_GPU_DEVICE (int box_no, int ix, int iy, int iz)
        {
            amrex::CellData<amrex::Real const> state = const_soln_arrays[box_no].cellData(i,j,k);
            amrex::CellData<amrex::Real> rhs = rhs_arrays[box_no].cellData(ix,iy,iz);
            ccz4rhs.compute(rhs, state);
        });
#endif
    }

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

    GRParmParse pp;
    bool activate_extraction{};
    pp.get("activate_extraction", activate_extraction);

    SphericalExtraction::params_t extraction_params; // TODO: Remove once extraction is fixed
    ExtractionTagger extraction_tagger(Geom().CellSize(0), Level(),
                                       extraction_params,
                                       activate_extraction);

    constexpr auto num_puncture_coords =
        static_cast<std::size_t>(AMREX_SPACEDIM * num_punctures);
    std::array<amrex::Real, num_puncture_coords> puncture_coords{};

    if (get_bhamr_ptr()->puncture_tracking_enabled)
    {
        puncture_coords = get_puncture_tracker().get_puncture_coords();
    }

    double bh1_mass;
    double bh2_mass;
    pp.get("bh1.mass", bh1_mass);
    pp.get("bh2.mass", bh2_mass);

    PunctureTagger<num_punctures> puncture_tagger(
        Geom().CellSize(0), Level(), get_gramr_ptr()->maxLevel(),
        puncture_coords, {bh1_mass, bh2_mass});

    amrex::ParallelFor(state_new, amrex::IntVect(0),
                       [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
                       {
                           chi_tagger(ix, iy, iz, tag_arrays[box_no],
                                      state_const_arrays[box_no]);

                           extraction_tagger(ix, iy, iz, tag_arrays[box_no]);

                           if (get_bhamr_ptr()->puncture_tracking_enabled)
                           {
                               puncture_tagger(ix, iy, iz, tag_arrays[box_no]);
                           }
                       });

    amrex::Gpu::streamSynchronize();
}

void BinaryBHLevel::specific_post_init()
{
    BL_PROFILE("BinaryBHLevel::specific_post_init()");

    if (get_bhamr_ptr()->puncture_tracking_enabled)
    {
        get_puncture_tracker().start_from_initial_punctures();
    }
}

void BinaryBHLevel::specific_post_restart()
{
    BL_PROFILE("BinaryBHLevel::specific_post_restart()");

    if (get_bhamr_ptr()->puncture_tracking_enabled)
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
    if (get_bhamr_ptr()->puncture_tracking_enabled)
    {
        get_puncture_tracker().write_plotfile(a_dir);
    }
}

void BinaryBHLevel::specific_post_checkpoint(const std::string &a_chk_dir,
                                             std::ostream & /*a_os*/)
{
    if (get_bhamr_ptr()->puncture_tracking_enabled)
    {
        get_puncture_tracker().checkpoint(a_chk_dir);
    }
}

void BinaryBHLevel::specificPostTimeStep()
{
    GRParmParse pp;
    int puncture_tracking_level;
    pp.get("puncture_tracking.level", puncture_tracking_level);
    int puncture_tracking_writeout_level;
    pp.get("puncture_tracking.writeout_level",
           puncture_tracking_writeout_level);

    // do puncture tracking on requested level
    if (get_bhamr_ptr()->puncture_tracking_enabled &&
        Level() == puncture_tracking_level)
    {
        BL_PROFILE("PunctureTracking");

        // only do the write out when we're at at a multiple of the
        // writeout_level
        bool write_punctures =
            at_level_timestep_multiple(puncture_tracking_writeout_level);
        amrex::Real current_time = get_state_data(state_index).curTime();
        amrex::Real dt           = get_gramr_ptr()->dtLevel(Level());
        get_puncture_tracker().track(current_time, dt, write_punctures);
    }
#if 0
//xxxxx specificPostTimeStep
    BL_PROFILE("BinaryBHLevel::specificPostTimeStep");

    bool first_step =
        (m_time == 0.); // this form is used when 'specificPostTimeStep' was
                        // called during setup at t=0 from Main
    // bool first_step = (m_time == m_dt); // if not called in Main

    bool activate_extraction{};
    pp.get("activate_extraction", activate_extraction);

    if (activate_extraction == 1)
    {
        int min_level = m_p.extraction_params.min_extraction_level();
        bool calculate_weyl = at_level_timestep_multiple(min_level);
        if (calculate_weyl)
        {
            // Populate the Weyl Scalar values on the grid
            fillAllGhosts();
            BoxLoops::loop(
                Weyl4(m_dx,),
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

    bool calculate_constraint_norms{};
    pp.get("calculate_constraint_norms", calculate_constraint_norms);


    if (calculate_constraint_norms)
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
