/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// GRTeclyn and test headers
#include "PunctureTrackerLevel.hpp"
#include "PunctureTagger.hpp"

// doctest header
#include "doctest.h"

BHAmr<PunctureTrackerLevel::num_punctures> *
PunctureTrackerLevel::get_bhamr_ptr()
{
    return dynamic_cast<BHAmr<num_punctures> *>(get_gr_amr_ptr());
}

PunctureTracker<PunctureTrackerLevel::num_punctures> &
PunctureTrackerLevel::get_puncture_tracker()
{
    return get_bhamr_ptr()->get_puncture_tracker();
}

void PunctureTrackerLevel::variableSetUp() { stateVariableSetUp(); }

void PunctureTrackerLevel::initData()
{
    amrex::MultiFab &state = get_new_data(state_index);
    const auto &arrs       = state.arrays();
    amrex::ParallelFor(state, state.nGrowVect(),
                       // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       {
                           const auto &array = arrs[box_no];
                           for (int icomp = 0; icomp < array.nComp(); ++icomp)
                           {
                               array(i, j, k, icomp) = 0.0;
                           }
                           array(i, j, k, c_shift2) = shift_y_val;
                       });
    amrex::Gpu::streamSynchronize();

    if (Level() == 0)
    {
        // need to set the puncture coordinates as we use it for the puncture
        // tagging
        GRParmParse puncture_tracking_pp("puncture_tracking");
        std::array<amrex::Real, AMREX_SPACEDIM * 2UL> initial_puncture_coords;
        puncture_tracking_pp.get("initial_coords", initial_puncture_coords);

        get_puncture_tracker().set_puncture_coords(initial_puncture_coords);
        // can't call start_from_initial_punctures() because we need the full
        // AMR grid first
    }
}

// Calculate RHS during RK4 substeps
// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
void PunctureTrackerLevel::specificEvalRHS(amrex::MultiFab &a_soln,
                                           amrex::MultiFab &a_rhs,
                                           const amrex::Real /*a_time*/)
{
    // We don't need any evolution in this test.
    a_rhs.setVal(0.0);
}

void PunctureTrackerLevel::tag_cells(amrex::TagBoxArray &a_tag_box_array,
                                     amrex::Real a_regrid_threshold)
{
    amrex::MultiFab &state_new = get_new_data(state_index);

    const auto &tag_arrs       = a_tag_box_array.arrays();
    const auto &state_new_arrs = state_new.const_arrays();

    std::array<amrex::Real, num_puncture_coords> puncture_coords{};

    puncture_coords = get_puncture_tracker().get_puncture_coords();

    amrex::Real fake_bh1_mass{};
    amrex::Real fake_bh2_mass{};

    GRParmParse test_pp("test");
    test_pp.get("fake_bh1_mass", fake_bh1_mass);
    test_pp.get("fake_bh2_mass", fake_bh2_mass);

    PunctureTagger<num_punctures> puncture_tagger(
        Geom().CellSize(0), Level(), get_gr_amr_ptr()->maxLevel(),
        puncture_coords, {fake_bh1_mass, fake_bh2_mass});

    amrex::ParallelFor(state_new, amrex::IntVect(0),
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       { puncture_tagger(i, j, k, tag_arrs[box_no]); });
    amrex::Gpu::streamSynchronize();
}

void PunctureTrackerLevel::specific_post_init()
{
    get_puncture_tracker().start_from_initial_punctures();
}

void PunctureTrackerLevel::specific_post_restart()
{

    std::string restart_checkpoint{};
    GRParmParse pp("amr");
    pp.get("restart", restart_checkpoint);
    get_puncture_tracker().restart(restart_checkpoint);
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
void PunctureTrackerLevel::specific_post_regrid(int a_lbase, int a_new_finest)
{
    // During initial data construction, we expect the finest level to increment
    // up to max_level so only do this for later steps
    if (get_gr_amr_ptr()->cumTime() > 0.0)
    {
        // Check the finest level is max_level
        CHECK(a_new_finest == get_gr_amr_ptr()->maxLevel());

        if (Level() > get_gr_amr_ptr()->maxLevel() - 2)
        {
            check_puncture_tagging();
        }
    }
}

void PunctureTrackerLevel::check_puncture_tagging()
{
    const int num_points_theta = 8;
    const int num_points_phi   = 8;

    amrex::Real fake_bh1_mass{};
    amrex::Real fake_bh2_mass{};

    GRParmParse test_pp("test");
    test_pp.get("fake_bh1_mass", fake_bh1_mass);
    test_pp.get("fake_bh2_mass", fake_bh2_mass);

    std::array<amrex::Real, num_punctures> fake_masses{fake_bh1_mass,
                                                       fake_bh2_mass};
    GRParmParse puncture_tagging_pp("puncture_tagging");
    amrex::Real level_separation{};
    puncture_tagging_pp.get("level_separation", level_separation);
    amrex::Real finest_level_factor{};
    puncture_tagging_pp.get("finest_level_factor", finest_level_factor);

    const int max_level        = get_gr_amr_ptr()->maxLevel();
    const amrex::Real exponent = max_level - Level();
    const amrex::Real factor =
        finest_level_factor * std::pow(level_separation, exponent);

    const auto &puncture_coords = get_puncture_tracker().get_puncture_coords();
    const auto &state_new       = get_new_data(state_index);
    const auto &box_array       = state_new.boxArray();

    for (int ipuncture = 0; ipuncture < num_punctures; ++ipuncture)
    {
        // test if sphere of points around the punctures are in the
        // BoxArray
        for (int itheta = 0; itheta < num_points_theta; ++itheta)
        {
            amrex::Real theta = (static_cast<amrex::Real>(itheta) + 0.5) *
                                M_PI / num_points_theta;
            for (int iphi = 0; iphi < num_points_phi; ++iphi)
            {
                amrex::Real phi = 2.0 * M_PI * static_cast<amrex::Real>(iphi) /
                                  num_points_phi;

                amrex::Real sphere_x = factor * fake_masses[ipuncture] *
                                           std::sin(theta) * std::cos(phi) +
                                       puncture_coords[ipuncture + 0];
                amrex::Real sphere_y = factor * fake_masses[ipuncture] *
                                           std::sin(theta) * std::sin(phi) +
                                       puncture_coords[ipuncture + 1];
                amrex::Real sphere_z =
                    factor * fake_masses[ipuncture] * std::cos(theta) +
                    puncture_coords[ipuncture + 2];

                amrex::RealVect sphere_coords{sphere_x, sphere_y, sphere_z};

                // skip points that are outside the problem domain - there will
                // be some if using reflective BCs
                if (!Geom().ProbDomain().contains(sphere_coords))
                {
                    continue;
                }

                bool point_covered_by_level = false;

                for (int ibox = 0; ibox < box_array.size(); ++ibox)
                {
                    const auto &box = box_array[ibox];
                    amrex::RealBox real_box{box, Geom().CellSize(),
                                            Geom().ProbLo()};
                    if (real_box.contains(sphere_coords))
                    {
                        point_covered_by_level = true;
                        break;
                    }
                }

                // This check might fail depending on the parameters if proper
                // nesting conditions prevents some cells from being refined
                INFO("Level: " << Level());
                INFO("Coordinates: (" << sphere_x << ", " << sphere_y << ", "
                                      << sphere_z << ")");
                INFO("Nearest puncture: ("
                     << puncture_coords[ipuncture + 0] << ", "
                     << puncture_coords[ipuncture + 1] << ", "
                     << puncture_coords[ipuncture + 2] << ")");
                CHECK(point_covered_by_level);
            }
        }
    }
}

void PunctureTrackerLevel::specific_post_checkpoint(
    const std::string &a_chk_dir, std::ostream & /*a_os*/)
{
    get_puncture_tracker().checkpoint(a_chk_dir);
}

void PunctureTrackerLevel::specificPostTimeStep()
{
    GRParmParse puncture_tracking_pp("puncture_tracking");
    int puncture_tracking_level{};
    puncture_tracking_pp.get("level", puncture_tracking_level);

    if (Level() == puncture_tracking_level)
    {
        bool write_punctures = false;
        amrex::Real cur_time = get_state_data(state_index).curTime();
        amrex::Real dt       = get_gr_amr_ptr()->dtLevel(Level());
        get_puncture_tracker().track(cur_time, dt, write_punctures);

        GRParmParse puncture_tracking_pp("puncture_tracking");
        std::array<amrex::Real, AMREX_SPACEDIM * 2UL> correct_puncture_coords;
        puncture_tracking_pp.get("initial_coords", correct_puncture_coords);

        for (int ipuncture = 0; ipuncture < num_punctures; ++ipuncture)
        {
            correct_puncture_coords[ipuncture * AMREX_SPACEDIM + 1] -=
                cur_time * shift_y_val;
        }
        auto puncture_coords = get_puncture_tracker().get_puncture_coords();
        constexpr int cout_precision = 16;
        for (int icoord = 0; icoord < num_puncture_coords; ++icoord)
        {
            INFO("puncture_coords["
                 << icoord << "] = " << std::setprecision(cout_precision)
                 << puncture_coords[icoord]);
            INFO("correct_puncture_coords["
                 << icoord << "] = " << std::setprecision(cout_precision)
                 << correct_puncture_coords[icoord]);
            CHECK(puncture_coords[icoord] ==
                  doctest::Approx(correct_puncture_coords[icoord])
                      .epsilon(1e-12));
        }
    }
}
