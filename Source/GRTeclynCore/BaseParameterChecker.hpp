/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef BASEPARAMETERCHECKER_HPP_
#define BASEPARAMETERCHECKER_HPP_

// General includes
#include "ArrayTools.hpp"
#include "BoundaryConditions.hpp"
#include "FilesystemTools.hpp"
#include "GRParmParse.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <unistd.h> // gives 'access'

class BaseParameterChecker
{
  public:
    BaseParameterChecker() = delete;

    // This function intentionally centralizes the validation of shared
    // parameters; splitting it would obscure dependencies between defaults.
    // NOLINTNEXTLINE(readability-function-cognitive-complexity)
    static void check_params()
    {
        GRParmParse amr_pp("amr");
        GRParmParse ccz4_pp("ccz4");
        GRParmParse evolution_pp("evolution");
        GRParmParse grteclyn_pp("grteclyn");
        GRParmParse particle_interpolator_pp("particle_interpolator");
        GRParmParse tagging_pp("tagging");

        // Grid setup

        int max_spatial_derivative_order = 4;
        evolution_pp.queryAdd("spatial_derivative_order",
                              max_spatial_derivative_order);
        if (max_spatial_derivative_order != 4 &&
            max_spatial_derivative_order != 6)
        {
            evolution_pp.error("spatial_derivative_order",
                               "only 4 and 6 are supported");
        }

        int max_grid_size = 64;
        amr_pp.queryAdd("max_grid_size", max_grid_size);

        if (max_grid_size < 0)
        {
            amr_pp.error("max_grid_size", "must be >= 0");
        }

        int blocking_factor = 8;
        amr_pp.queryAdd("blocking_factor", blocking_factor);

        if (blocking_factor < 1)
        {
            amr_pp.error("blocking_factor", "must be >= 1");
        }

        if (max_grid_size % blocking_factor != 0)
        {
            amr_pp.error("blocking_factor",
                         "must divide max_grid_size/max_box_size");
        }

        int num_ghosts = (max_spatial_derivative_order == 6) ? 4 : 3;
        evolution_pp.queryAdd("num_ghosts", num_ghosts);

        // the following check assumes you will be taking one-sided derivatives
        // of the order given by max_spatial_derivative_order

        if ((num_ghosts < ((max_spatial_derivative_order == 6) ? 4 : 3)) ||
            (num_ghosts > blocking_factor))
        {
            evolution_pp.error(
                "num_ghosts",
                "must be >= 3 (4th order derivatives) or 4 (6th order "
                "derivatives) and <= blocking_factor");
        }

        // check the restart_file exists and can be read if restarting from a
        // checkpoint
        if (amr_pp.contains("restart"))
        {
            std::string restart_file;
            amr_pp.get("restart", restart_file);
            if (access((restart_file).c_str(), R_OK) != 0)
            {
                amr_pp.error("restart", "file cannot be opened for reading");
            }
        }

        int n_error_buf = 3; // Amount the tagged region is grown by

        amr_pp.queryAdd("n_error_buf", n_error_buf);
        if (n_error_buf < 0)
        {
            amr_pp.error("n_error_buf", "must be >= 0");
        }

        // Multiplier for the proper-nesting buffer. AMReX measures this in
        // coarse-cell blocks of max(1, blocking_factor / ref_ratio).
        int n_proper = 1;
        amr_pp.queryAdd("n_proper", n_proper);
        // GRTeclyn fixes the refinement ratio to 2 below.
        constexpr int fixed_ref_ratio = 2;
        const int proper_nesting_buffer =
            n_proper * std::max(1, blocking_factor / fixed_ref_ratio);
        if (proper_nesting_buffer < num_ghosts)
        {
            amr_pp.error("n_proper",
                         "times max(1, blocking_factor / ref_ratio) must be >= "
                         "num_ghosts for proper nesting");
        }

        // Minimum fraction of tagged cells required in each generated grid.
        double grid_eff = 0.7;
        amr_pp.queryAdd("grid_eff", grid_eff);
        if (grid_eff <= 0.0 || grid_eff > 1.0)
        {
            amr_pp.error("grid_eff", "must be > 0 and <= 1");
        }

        double dt_multiplier = 0.25;
        evolution_pp.queryAdd("dt_multiplier", dt_multiplier);
        if (dt_multiplier <= 0.0)
        {
            evolution_pp.error("dt_multiplier", "must be > 0.0");
        }
        else if (dt_multiplier >= 1.0)
        {
            evolution_pp.error("dt_multiplier", "must be < 1.0 for stability");
        }
        else if (dt_multiplier > 0.5)
        {
            evolution_pp.warning("dt_multiplier",
                                 "is unlikely to be stable for > 0.5");
        }

        // For n_cell and prob_extent, must factor in reflective boundaries
        std::array<int, AMREX_SPACEDIM> n_cell;
        amr_pp.get("n_cell", n_cell);

        GRParmParse geom_pp("geometry");

        std::array<double, AMREX_SPACEDIM> prob_extent{};
        geom_pp.get("prob_extent", prob_extent);

        std::array<double, AMREX_SPACEDIM> dx{};
        double dx_tol = 1e-10;
        FOR (idir)
        {
            if (prob_extent[idir] <= 0.0)
            {
                geom_pp.error("prob_extent", "must be > 0 in all directions");
            }
            dx[idir] = prob_extent[idir] / n_cell[idir];
        }
        FOR (idir)
        {
            if (std::abs(dx[idir] - dx[(idir + 2) % 3]) > dx_tol)
            {
                geom_pp.error("prob_extent",
                              "does not give equal dx in each direction with "
                              "provided amr.n_cell");
            }
        }
        double coarsest_dx = prob_extent[0] / n_cell[0];

        geom_pp.add("coarsest_dx", coarsest_dx);

        std::array<int, AMREX_SPACEDIM> is_periodic = {0, 0, 0};
        geom_pp.queryAdd("is_periodic", is_periodic);

        // Periodicity and boundaries
        // TODO: Improve so params don't need to be filled
        BoundaryConditions::params_t::check_params();
        BoundaryConditions::params_t boundary_params;
        boundary_params.fill_params();

        // Work out the default center, factoring in reflective boundaries

        for (const char *prob_bound_name : {"prob_lo", "prob_hi"})
        {
            if (geom_pp.contains(prob_bound_name))
            {
                geom_pp.warning(
                    prob_bound_name,
                    "not implemented, assumed to be (0,0,0) and prob_extent");
            }
        }

        std::array<double, AMREX_SPACEDIM> center{};
        FOR (idir)
        {
            if ((boundary_params.lo_condition[idir] ==
                 BoundaryConditions::REFLECTIVE_BC) &&
                (boundary_params.hi_condition[idir] !=
                 BoundaryConditions::REFLECTIVE_BC))
            {
                center[idir] = 0.;
            }
            else if ((boundary_params.hi_condition[idir] ==
                      BoundaryConditions::REFLECTIVE_BC) &&
                     (boundary_params.lo_condition[idir] !=
                      BoundaryConditions::REFLECTIVE_BC))
            {
                center[idir] = prob_extent[idir];
            }
            else
            {
                center[idir] = 0.5 * prob_extent[idir];
            }
        }
        geom_pp.queryAdd("center", center);

        int max_level = 0; // the max number of regriddings to do
        amr_pp.queryAdd("max_level", max_level);
        if (max_level < 0)
        {
            amr_pp.error("max_level", "must be >= 0");
        }

        // The refinement ratio is fixed to 2 since GRTeclyn is not tested with
        // other values.
        amrex::Vector<int> ref_ratio(max_level, 2); // ref ratios between levels
        amr_pp.addarr("ref_ratio", ref_ratio);

        // Regridding interval on each level, with size max_level (i.e.
        // num_levels-1) since regridding on max level does nothing
        amrex::Vector<int> regrid_int(max_level, -1);
        if (max_level > 0)
        {
            regrid_int[0] = 1;
        }
        amr_pp.queryAdd("regrid_int", regrid_int);

        if (!tagging_pp.contains("thresholds"))
        {
            double regrid_threshold = 0.5;
            tagging_pp.queryAdd("threshold", regrid_threshold);

            amrex::Vector<double> regrid_thresholds(max_level + 1,
                                                    regrid_threshold);
            tagging_pp.queryAdd("thresholds", regrid_thresholds);
        }

        int check_int = -1; // Steps between checkpoint file outputs
        amr_pp.queryAdd("check_int", check_int);

        int plot_int = -1; // Steps between plot file outputs
        amr_pp.queryAdd("plot_int", plot_int);

        double stop_time = 1.; // The stop time
        evolution_pp.queryAdd("stop_time", stop_time);
        if (stop_time < 0.0)
        {
            evolution_pp.warning("stop_time",
                                 "is negative, disabling this stopping "
                                 "criterion");
        }

        int max_steps = 1000000;
        evolution_pp.queryAdd("max_steps", max_steps);
        if (max_steps < 0)
        {
            evolution_pp.warning("max_steps",
                                 "is negative, disabling this stopping "
                                 "criterion");
        }

        FOR (idir)
        {
            if (n_cell[idir] % blocking_factor != 0)
            {
                amr_pp.error("n_cell", "must divide blocking_factor");
            }
        }

        double sigma = 0.1; // Kreiss-Oliger dissipation parameter
        // Dissipation
        evolution_pp.queryAdd("sigma", sigma);
        if (sigma < 0.0 || sigma > 2.0 / dt_multiplier)
        {
            evolution_pp.warning("sigma",
                                 "must be >= 0.0 and <= 2 / dt_multiplier for "
                                 "stability (see Alcubierre p344)");
        }

        // Nan Check and min chi and lapse values
        bool nan_check = true;
        evolution_pp.queryAdd("nan_check", nan_check);
        if (!nan_check)
        {
            evolution_pp.warning("nan_check",
                                 "should not normally be disabled");
        }

        double min_chi   = 1e-4;
        double min_lapse = 1e-4;
        ccz4_pp.queryAdd("min_chi", min_chi);
        ccz4_pp.queryAdd("min_lapse", min_lapse);
        if (min_chi < 0.0)
        {
            ccz4_pp.warning("min_chi", "should normally be >= 0.0");
        }
        if (min_lapse < 0.0)
        {
            ccz4_pp.warning("min_lapse", "should normally be >= 0.0");
        }

        bool particle_interpolator_verbosity = false;
        particle_interpolator_pp.queryAdd("verbosity",
                                          particle_interpolator_verbosity);

        std::string output_path = ".";
        grteclyn_pp.queryAdd("output_path", output_path);

        std::string plot_directory = output_path + "/plots";
        if (!FilesystemTools::directory_exists(plot_directory))
        {
            FilesystemTools::mkdir_recursive(plot_directory);
        }

        std::string checkpoint_directory = output_path + "/checkpoints";
        if (!FilesystemTools::directory_exists(checkpoint_directory))
        {
            FilesystemTools::mkdir_recursive(checkpoint_directory);
        }

        std::string plot_file  = plot_directory + "/plt";
        std::string check_file = checkpoint_directory + "/chk";
        amr_pp.queryAdd("plot_file", plot_file);
        amr_pp.queryAdd("check_file", check_file);
    }
};

#endif /* BASEPARAMETERCHECKER_HPP_ */
