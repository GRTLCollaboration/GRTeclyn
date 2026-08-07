/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef BASEPARAMETERCHECKER_HPP_
#define BASEPARAMETERCHECKER_HPP_

// General includes
#include "ArrayTools.hpp"
#include "BoundaryConditions.hpp"
#include "CCZ4RHS.hpp"
#include "FilesystemTools.hpp"
#include "GRParmParse.hpp"
#include "SphericalExtraction.hpp"
#include "StateVariables.hpp"
#include "VariableType.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <unistd.h> // gives 'access'

// add this type alias here for backwards compatibility
using extraction_params_t = SphericalExtraction::params_t;

class BaseParameterChecker
{
  public:
    BaseParameterChecker()
    {
        check_amrex_params();

        check_grteclyn_params();
    }

    void check_amrex_params()
    {
        GRParmParse amr_pp("amr");

        // Grid setup

        int max_spatial_derivative_order = 4;
        amr_pp.queryAdd("max_spatial_derivative_order",
                        max_spatial_derivative_order);
        if (!(max_spatial_derivative_order == 4 ||
              max_spatial_derivative_order == 6))
        {
            amr_pp.error("max_spatial_derivative_order",
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
        amr_pp.queryAdd("num_ghosts", num_ghosts);

        // the following check assumes you will be taking one-sided derivatives
        // of the order given by max_spatial_derivative_order

        if ((num_ghosts < ((max_spatial_derivative_order == 6) ? 4 : 3)) ||
            (num_ghosts > blocking_factor))
        {
            amr_pp.error("num_ghosts",
                         "must be >= 3 (4th order derivatives) or 4 (6th order "
                         "derivatives) and <= blocking_factor");
        }

        // check the restart_file exists and can be read if restarting from a
        // checkpoint
        if (amr_pp.contains("restart"))
        {
            std::string restart_file;
            amr_pp.load("restart", restart_file);
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

        int n_proper = 1;
        amr_pp.queryAdd("n_proper", n_proper);
        // assume ref_ratio is always 2
        if (blocking_factor * n_proper < num_ghosts)
        {
            amr_pp.error("n_proper", " times blocking_factor must be >= "
                                     "num_ghosts for proper nesting");
        }

        double grid_eff =
            0.7; // determines how fussy the regridding is about tags
        amr_pp.queryAdd("grid_eff", grid_eff);
        if (grid_eff <= 0.0 || grid_eff > 1.0)
        {
            amr_pp.error("grid_eff", "must be > 0 and <= 1");
        }

        double dt_multiplier = 0.25;
        amr_pp.queryAdd("dt_multiplier", dt_multiplier);
        if (dt_multiplier <= 0.0)
        {
            amr_pp.error("dt_multiplier", "must be > 0.0");
        }
        else if (dt_multiplier >= 1.0)
        {
            amr_pp.error("dt_multiplier", "must be < 1.0 for stability");
        }
        else if (dt_multiplier > 0.5)
        {
            amr_pp.warning("dt_multiplier",
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

        amr_pp.add("coarsest_dx", coarsest_dx);

        std::array<int, AMREX_SPACEDIM> is_periodic = {0, 0, 0};
        geom_pp.queryAdd("is_periodic", is_periodic);

        // Periodicity and boundaries
        BoundaryConditions::params_t::check_params();

        boundary_params.fill_params();

        // Work out the default center, factoring in reflective boundaries
        // Assumes prob_lo = 0 0 0 (the default)

        std::array<double, AMREX_SPACEDIM> center{};
        FOR (idir)
        {
            if ((boundary_params.lo_boundary[idir] ==
                 BoundaryConditions::REFLECTIVE_BC) &&
                (boundary_params.hi_boundary[idir] !=
                 BoundaryConditions::REFLECTIVE_BC))
            {
                center[idir] = 0.;
            }
            else if ((boundary_params.hi_boundary[idir] ==
                      BoundaryConditions::REFLECTIVE_BC) &&
                     (boundary_params.lo_boundary[idir] !=
                      BoundaryConditions::REFLECTIVE_BC))
            {
                center[idir] = prob_extent[idir];
            }
            else
            {
                center[idir] = 0.5 * prob_extent[idir];
            }
        }
        amr_pp.queryAdd("center", center);

        int max_level = 0; // the max number of regriddings to do
        amr_pp.queryAdd("max_level", max_level);
        if (max_level < 0)
        {
            amr_pp.error("max_level", "must be >= 0");
        }

        // the reference ratio is hard coded to 2
        // in principle it can be set to other values, but this is
        // not recommended since we do not test GRChombo with other
        // refinement ratios - use other values at your own risk
        amrex::Vector<int> ref_ratio(max_level, 2); // ref ratios between levels
        amr_pp.queryAdd("ref_ratio", ref_ratio);

        // Regridding interval on each level, with size max_level (i.e.
        // num_levels-1) since regridding on max level does nothing
        amrex::Vector<int> regrid_int(max_level,
                                      2); // steps between regrid at each level
        amr_pp.queryAdd("regrid_int", regrid_int);

        if (amr_pp.contains("regrid_thresholds"))
        {
            amrex::Print() << "Using multiple regrid thresholds." << '\n';
        }
        else
        {
            amrex::Print() << "Using single regrid threshold." << '\n';

            double regrid_threshold = 0.5;
            amr_pp.queryAdd("regrid_threshold", regrid_threshold);

            amrex::Vector<double> regrid_thresholds(max_level + 1,
                                                    regrid_threshold);
            amr_pp.queryAdd("regrid_thresholds", regrid_thresholds);
        }

        int check_int = 1; // Steps between checkpoint file outputs
        amr_pp.queryAdd("check_int", check_int);

        int plot_int = 0; // Steps between plot file outputs
        amr_pp.queryAdd("plot_int", plot_int);

        double stop_time = 1.; // The stop time
        amr_pp.queryAdd("stop_time", stop_time);

        int max_steps = 1000000;
        amr_pp.queryAdd("max_steps", max_steps);

        FOR (idir)
        {
            if (n_cell[idir] % blocking_factor != 0)
            {
                amr_pp.error("n_cell", "must divide blocking_factor");
            }
        }

        // TODO: Check extraction (and interpolation) params e.g. origin and
        // reflective_domain_lo/hi
    }

    void check_grteclyn_params()
    {
        GRParmParse pp;
        double dt_multiplier;
        pp.get("amr.dt_multiplier", dt_multiplier);

        double sigma = 0.1; // Kreiss-Oliger dissipation parameter
        // Dissipation
        pp.queryAdd("sigma", sigma);
        if (sigma < 0.0 || sigma > 2.0 / dt_multiplier)
        {
            pp.warning("sigma", "must be >= 0.0 and <= 2 / dt_multiplier for "
                                "stability (see Alcubierre p344)");
        }

        // Nan Check and min chi and lapse values
        bool nan_check = true;
        pp.queryAdd("nan_check", nan_check);
        if (!nan_check)
        {
            pp.warning("nan_check", "should not normally be disabled");
        }

        double min_chi   = 1e-4;
        double min_lapse = 1e-4;
        pp.queryAdd("min_chi", min_chi);
        pp.queryAdd("min_lapse", min_lapse);

        // not sure these are necessary hence commented out
        // check_parameter("min_chi", min_chi, (min_chi >= 0.0), "must be >=
        // 0.0"); check_parameter("min_lapse", min_lapse, (min_lapse >= 0.0)
        // "must be >= 0.0");

        int formulation = CCZ4RHS<>::USE_CCZ4; // Whether to use BSSN or CCZ4
        pp.queryAdd("formulation", formulation);

        if (formulation != CCZ4RHS<>::USE_CCZ4 &&
            formulation != CCZ4RHS<>::USE_BSSN)
        {
            pp.error("formulation", "must be 0 or 1");
        }

        if (formulation == CCZ4RHS<>::USE_CCZ4)
        {
            CCZ4_params_t::check_params();
        }
        else if (formulation == CCZ4RHS<>::USE_BSSN)
        {
            if (pp.contains("ccz4.kappa1") || pp.contains("ccz4.kappa2") ||
                pp.contains("ccz4.kappa3"))
            {
                pp.warning("kappa1/2/3",
                           "should not be provided with BSSN formulation, "
                           "setting them all to zero");
            }
            pp.add("ccz4.kappa1", 0.0);
            pp.add("ccz4.kappa2", 0.0);
            pp.add("ccz4.kappa3", 0.0);
        }

        // Extraction params
        bool activate_extraction = false;
        pp.queryAdd("activate_extraction", activate_extraction);

        if (activate_extraction)
        {
            SphericalExtraction::params_t::check_params();
        }

        // TODO: Remove once extraction is fixed
        std::array<double, AMREX_SPACEDIM> center{};
        pp.get("amr.center", center);
        std::array<double, AMREX_SPACEDIM> extraction_center = center;
        pp.queryAdd("extraction_center", extraction_center);
    }

    // General parameters

    // boundaries.

    // Boundary conditions
    BoundaryConditions::params_t boundary_params; // set boundaries in each dir

  protected:
    // the low and high corners of the domain taking into account reflective BCs
    // only used in parameter checks hence protected
    std::array<double, AMREX_SPACEDIM> reflective_domain_lo{},
        reflective_domain_hi{};

    // use this error function instead of MayDay::error as this will only
    // print from rank 0
    static void error(const std::string &a_error_message)
    {
        if (amrex::ParallelDescriptor::MyProc() == 0)
        {
            amrex::Abort(a_error_message.c_str());
        }
    }

    template <typename T>
    void check_parameter(const std::string &a_name, const T &a_value,
                         const bool a_valid,
                         const std::string &a_invalid_explanation)
    {
        if (a_valid)
        {
            return;
        }
        std::ostringstream error_message_ss;
        error_message_ss << "Parameter: " << a_name << " = " << a_value
                         << " is invalid: " << a_invalid_explanation;
        error(error_message_ss.str());
    }

    template <typename T>
    void warn_parameter(const std::string &a_name, T a_value,
                        const bool a_nowarn,
                        const std::string &a_warning_explanation)
    {
        if (a_nowarn)
        {
            return;
        }
        // only print the warning from rank 0
        if (amrex::ParallelDescriptor::MyProc() == 0)
        {
            std::ostringstream warning_message_ss;
            warning_message_ss << "Parameter: " << a_name << " = " << a_value
                               << " warning: " << a_warning_explanation;
            amrex::Warning(warning_message_ss.str().c_str());
        }
    }

    template <typename T, size_t N>
    void check_array_parameter(const std::string &a_name,
                               const std::array<T, N> &a_value,
                               const bool a_valid,
                               const std::string &a_invalid_explanation)
    {
        std::string value_str = ArrayTools::to_string(a_value);
        check_parameter(a_name, value_str, a_valid, a_invalid_explanation);
    }

    template <typename T, size_t N>
    void warn_array_parameter(const std::string &a_name,
                              const std::array<T, N> &a_value,
                              const bool a_nowarn,
                              const std::string &a_warning_explanation)
    {
        std::string value_str = ArrayTools::to_string(a_value);
        check_parameter(a_name, value_str, a_nowarn, a_warning_explanation);
    }
};

#endif /* CHOMBOPARAMETERS_HPP_ */
