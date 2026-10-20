/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SPECTRALWISHESAMRTESTLEVEL_HPP_
#define SPECTRALWISHESAMRTESTLEVEL_HPP_

#include "doctest.h"

#include "DefaultLevelBld.hpp"
#include "DefaultPotential.hpp"
#include "FixedGridsTagger.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRAmrLevel.hpp"
#include "ScalarField.hpp"
#include "SpectralWishes.hpp"

#include <filesystem>
#include <string>

/*!
 * Two levels (the finer one placed by FixedGridsTagger) with analytic data,
 * phi = 1, Pi = x and chi = 1/4, which is not evolved, to check SpectralWishes
 * over several levels against exact results. Each set of wishes in
 * params_amr_test.txt has its own parameter prefix.
 */
class SpectralWishesAMRTestLevel : public GRAmrLevel
{
  public:
    using GRAmrLevel::GRAmrLevel;

    using DefaultScalarField =
        ScalarField<DefaultPotential, FourthOrderDerivatives>;

    static constexpr amrex::Real s_phi = 1.0;
    static constexpr amrex::Real s_chi = 0.25;

    static void variableSetUp() { state_variable_set_up(); };

    void specific_advance() override {};

    void initData() override
    {
        const auto dx          = Geom().CellSizeArray();
        const auto prob_lo     = Geom().ProbLoArray();
        amrex::MultiFab &state = get_new_data(state_index);
        state.setVal(0.0);
        const auto &state_arrays = state.arrays();
        // NOLINTBEGIN(bugprone-easily-swappable-parameters)
        amrex::ParallelFor(
            state, state.nGrowVect(),
            [=] AMREX_GPU_DEVICE(int ibox, int ix, int iy, int iz)
            // NOLINTEND(bugprone-easily-swappable-parameters)
            {
                const auto &cell        = state_arrays[ibox];
                cell(ix, iy, iz, c_phi) = s_phi;
                cell(ix, iy, iz, c_Pi)  = prob_lo[0] + (ix + 0.5) * dx[0];
                cell(ix, iy, iz, c_chi) = s_chi;
            });
        amrex::Gpu::streamSynchronize();
    };

    // The data is not evolved
    void specific_eval_rhs(amrex::MultiFab & /*a_soln*/, amrex::MultiFab &a_rhs,
                           amrex::Real /*a_time*/) override
    {
        a_rhs.setVal(0.0);
    };

    void specific_update_ode(amrex::MultiFab & /*a_soln*/) override {};

    void tag_cells(amrex::TagBoxArray &a_tag_box_array,
                   amrex::Real /*a_regrid_threshold*/) final
    {
        const auto &tag_arrays = a_tag_box_array.arrays();
        const FixedGridsTagger tagger(Geom().CellSize(0), Level());
        amrex::ParallelFor(
            a_tag_box_array,
            [=] AMREX_GPU_DEVICE(int ibox, int ix, int iy, int iz)
            { tagger(ix, iy, iz, tag_arrays[ibox]); });
        amrex::Gpu::streamSynchronize();
    };

    void specific_post_timestep() override
    {
        if (Level() != 0)
        {
            return;
        }

        auto *gr_amr_ptr = get_gr_amr_ptr();
        REQUIRE(gr_amr_ptr->finestLevel() == 1);

        // Pi = x on the coarse grid of 32 cells over [0, L] and the fine
        // level, wherever FixedGridsTagger put it
        const amrex::Real length    = Geom().ProbLength(0);
        const amrex::Real coarse_dx = Geom().CellSize(0);
        const amrex::Real fine_dx   = gr_amr_ptr->Geom(1).CellSize(0);
        const amrex::Real volume    = Geom().ProbSize();

        // Moments of x over the fine level's cells (exact sums over its
        // BoxArray)
        amrex::Real fine_num_cells{0.};
        amrex::Real fine_sum_x{0.};
        amrex::Real fine_sum_x2{0.};
        const amrex::BoxArray &fine_box_array = gr_amr_ptr->boxArray(1);
        for (int ibox = 0; ibox < fine_box_array.size(); ++ibox)
        {
            const amrex::Box &box = fine_box_array[ibox];
            const auto cells_per_x =
                static_cast<amrex::Real>(box.numPts() / box.length(0));
            for (int ix = box.smallEnd(0); ix <= box.bigEnd(0); ++ix)
            {
                const amrex::Real x  = (ix + 0.5) * fine_dx;
                fine_num_cells      += cells_per_x;
                fine_sum_x          += cells_per_x * x;
                fine_sum_x2         += cells_per_x * x * x;
            }
        }
        const amrex::Real fine_volume = fine_num_cells * std::pow(fine_dx, 3);
        const amrex::Box fine_box     = fine_box_array.minimalBox();
        const amrex::Real fine_min_x  = (fine_box.smallEnd(0) + 0.5) * fine_dx;
        const amrex::Real fine_max_x  = (fine_box.bigEnd(0) + 0.5) * fine_dx;
        const amrex::Real fine_mean   = fine_sum_x / fine_num_cells;
        const amrex::Real fine_variance =
            fine_sum_x2 / fine_num_cells - fine_mean * fine_mean;

        // The refined region is not centred so counting covered cells twice
        // would move the composite mean away from L/2
        REQUIRE(fine_mean != doctest::Approx(0.5 * length));

        // The coarse cells covered by the fine level hold data interpolated
        // from it (in GRAmrLevel::post_timestep) which is not exactly x, so
        // compare the coarse level on its own with AMReX's sums of its data
        const amrex::MultiFab &coarse_state = get_new_data(state_index);
        const auto coarse_num_cells =
            static_cast<amrex::Real>(Geom().Domain().numPts());
        const amrex::Real coarse_mean =
            coarse_state.sum(c_Pi) / coarse_num_cells;
        const amrex::Real coarse_variance =
            std::pow(coarse_state.norm2(c_Pi), 2) / coarse_num_cells -
            coarse_mean * coarse_mean;

        // Over cells of width h the midpoint rule is exact for x and gives
        // int (x - L/2)^2 dV - V h^2 / 12 for the second moment, so on the
        // composite hierarchy
        const amrex::Real mean = 0.5 * length;
        const amrex::Real composite_variance =
            length * length / 12.0 -
            ((volume - fine_volume) * coarse_dx * coarse_dx +
             fine_volume * fine_dx * fine_dx) /
                (12.0 * volume);

        // Mean Variance of Pi on the calling level (default)
        check_wishes("wishes_default", {coarse_mean, coarse_variance});

        // Over both levels: Mean Variance VolumeMean of Pi, Integral
        // VolumeIntegral of phi (= 1, with sqrt(gamma) = chi^-3/2 = 8), Min Max
        // L1Norm of Pi, L2Norm of phi and Max of the magnitude phi:Pi. The
        // extreme values of x are on the coarse level next to the boundary.
        const amrex::Real sqrt_gamma = std::pow(s_chi, -1.5);
        const amrex::Real max_x      = length - 0.5 * coarse_dx;
        check_wishes("wishes_composite",
                     {mean, composite_variance, mean, s_phi * volume,
                      s_phi * sqrt_gamma * volume, 0.5 * coarse_dx, max_x, mean,
                      s_phi, std::sqrt(s_phi * s_phi + max_x * max_x)});

        // Mean Variance of Pi on each level
        check_wishes("wishes_each",
                     {coarse_mean, coarse_variance, fine_mean, fine_variance},
                     {"Mean_Pi_lev0", "Variance_Pi_lev0", "Mean_Pi_lev1",
                      "Variance_Pi_lev1"});

        // Mean Variance of Pi, Integral of phi and Min Max of Pi on the fine
        // level only
        check_wishes("wishes_fine",
                     {fine_mean, fine_variance, s_phi * fine_volume, fine_min_x,
                      fine_max_x});

        // Leaving out a layer of width 1 (2 coarse cells) next to every face,
        // which only covers coarse cells
        const amrex::Real excluded_width  = 2 * coarse_dx;
        const amrex::Real included_length = length - 2 * excluded_width;
        check_wishes("wishes_excluded",
                     {s_phi * std::pow(included_length, 3), mean,
                      excluded_width + 0.5 * coarse_dx,
                      length - excluded_width - 0.5 * coarse_dx});

        // On the fine level, leaving out x < 2.5 (5 coarse cells, so 10 fine
        // cells), which would not remove anything if the width were not
        // converted to fine cells
        const amrex::Real fine_excluded_x = 5 * coarse_dx;
        REQUIRE(fine_excluded_x > fine_min_x);
        check_wishes("wishes_fine_excluded",
                     {s_phi * fine_volume *
                          (fine_max_x + 0.5 * fine_dx - fine_excluded_x) /
                          (fine_max_x - fine_min_x + fine_dx),
                      fine_excluded_x + 0.5 * fine_dx});

        // Groups inherit max_level = -1 from wishes_groups (otherwise the
        // mean would be the coarse one, not exactly L/2) and interior adds
        // its own exclusion
        check_wishes("wishes_groups", {s_phi * volume, mean}, {}, "whole");
        check_wishes("wishes_groups",
                     {s_phi * std::pow(included_length, 3), mean}, {},
                     "interior");

        // write_groups writes one file per group (the path ./ is inherited)
        SpectralWishes<>::write_groups(
            state_index, get_state_data(state_index).curTime(), Level(),
            get_gr_amr_ptr(), "wishes_groups");
        if (amrex::ParallelDescriptor::IOProcessor())
        {
            for (const std::string group : {"whole", "interior"})
            {
                const std::string filename =
                    "./wishes_groups_" + group + ".dat";
                INFO("Group file ", filename);
                CHECK(std::filesystem::exists(filename));
                std::filesystem::remove(filename);
            }
        }
    };

  protected:
    //! Compute the wishes with prefix a_prefix (and group a_group if not
    //! empty) and compare them (and optionally their labels) with the
    //! expected values
    void check_wishes(const std::string &a_prefix,
                      const amrex::Vector<amrex::Real> &a_expected,
                      const amrex::Vector<std::string> &a_expected_labels = {},
                      const std::string &a_group                          = "")
    {
        const amrex::Real time = get_state_data(state_index).curTime();
        SpectralWishes<DefaultScalarField> wishes(state_index, time, Level(), 0,
                                                  a_prefix, a_group);
        const auto results = wishes.compute(get_gr_amr_ptr());
        const auto labels  = wishes.get_labels();

        INFO("Spectral wishes with prefix ", a_prefix, " and group ", a_group);
        if (!a_group.empty())
        {
            // One file per group, named after it
            CHECK(wishes.get_file_prefix().ends_with(a_prefix + "_" + a_group));
        }
        REQUIRE(results.size() == a_expected.size());
        if (!a_expected_labels.empty())
        {
            CHECK(labels == a_expected_labels);
        }
        for (amrex::Long iresult = 0; iresult < results.size(); ++iresult)
        {
            amrex::Print() << a_prefix << " " << labels[iresult] << ": "
                           << results[iresult] << " (expected "
                           << a_expected[iresult] << ")\n";
            INFO("Spectral wish: ", labels[iresult]);
            CHECK(results[iresult] ==
                  doctest::Approx(a_expected[iresult]).epsilon(1e-12));
        }
    }
};

#endif /* SPECTRALWISHESAMRTESTLEVEL_HPP_ */
