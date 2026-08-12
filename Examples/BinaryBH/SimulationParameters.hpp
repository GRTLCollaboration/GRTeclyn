/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP
#define SIMULATIONPARAMETERS_HPP

// General includes
#include "BaseParameterChecker.hpp"
#include "GRParmParse.hpp"

// Problem specific includes:
#include "ArrayTools.hpp"
#include "BoostedBHInitialData.hpp"
#include "PunctureTracker.hpp"
#include "CCZ4RHS.hpp"
#ifdef USE_TWOPUNCTURES
#include "TwoPuncturesInitialData.hpp"
#endif

class SimulationParameters
{
  public:
    // NOLINTNEXTLINE(readability-identifier-length)
    SimulationParameters() = delete;

    static void check_params()
    {
        BaseParameterChecker::check_params();
        
        read_shared_params();

#ifndef USE_TWOPUNCTURES
        BoostedBHInitialData::params_t::check_params(1);
        BoostedBHInitialData::params_t::check_params(2);
#endif
    }

    /// Read shared parameters
    // NOLINTNEXTLINE(readability-identifier-length)
    static void read_shared_params()
    {
        GRParmParse pp;
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
        
        // Do we want puncture tracking and constraint norm calculation?
        bool puncture_tracking_enabled{false};
        pp.queryAdd("puncture_tracking.enabled", puncture_tracking_enabled);
        if (puncture_tracking_enabled)
        {
            puncture_tracker_params_t::check_params();
        }

        bool calculate_constraint_norms = false;
        pp.queryAdd("calculate_constraint_norms", calculate_constraint_norms);

        int max_level;
        pp.get("amr.max_level", max_level);
        int puncture_tracking_level = max_level;
        pp.queryAdd("puncture_tracking_level", puncture_tracking_level);
        if (puncture_tracking_level < 0) || (puncture_tracking_level > max_level)
        {
            pp.error("puncture_tracking_level", "must be between 0 and max_level (inclusive)");
        }
    }

#if 0
    void check_tp_params()
    {
        // These checks are mostly taken from the Einstein Toolkit thorn
        // documentation:
        // https://einsteintoolkit.org/thornguide/EinsteinInitialData/TwoPunctures/documentation.html
        std::string mass_plus_name, mass_minus_name;
        if (tp_params.give_bare_mass)
        {
            mass_minus_name = "TP_mass_minus";
            mass_plus_name  = "TP_mass_plus";
        }
        else
        {
            mass_minus_name = "TP_target_mass_minus";
            mass_plus_name  = "TP_target_mass_plus";
            check_parameter("TP_adm_tol", tp_params.adm_tol,
                            tp_params.adm_tol > 0., "must be > 0.0");
        }
        check_parameter(mass_minus_name, bh1_params.mass, bh1_params.mass >= 0.,
                        "mustd be >= 0.0");
        check_parameter(mass_plus_name, bh2_params.mass, bh2_params.mass >= 0.,
                        "must be >= 0.0");

        int offset_dir = (!tp_params.swap_xz) ? 0 : 2;
        std::array<double, AMREX_SPACEDIM> center{};
        pp.get("amr.center", center);
        warn_parameter("TP_offset_minus", tp_offset_minus,
                       tp_offset_minus < (ivN[offset_dir] + 1) * coarsest_dx -
                                             center[offset_dir],
                       "should be within the computational domain");
        warn_parameter("TP_offset_plus", tp_offset_plus,
                       tp_offset_plus < (ivN[offset_dir] + 1) * coarsest_dx -
                                            center[offset_dir],
                       "should be within the computational domain");
        check_parameter("TP_npoints_A", tp_params.npoints_A,
                        tp_params.npoints_A >= 4, "must be >= 4");
        check_parameter("TP_npoints_B", tp_params.npoints_B,
                        tp_params.npoints_B >= 4, "must be >= 4");
        check_parameter("TP_npoints_phi", tp_params.npoints_phi,
                        tp_params.npoints_phi >= 4 &&
                            tp_params.npoints_phi % 2 == 0,
                        "must be >= 4 and divisible by 2");
        check_parameter("TP_Newton_maxit", tp_params.Newton_maxit,
                        tp_params.Newton_maxit >= 0, "must be >= 0");
        check_parameter("TP_Newton_tol", tp_params.Newton_tol,
                        tp_params.Newton_tol >= 0., "must be >= 0.0");
        check_parameter("TP_epsilon", tp_params.TP_epsilon,
                        tp_params.TP_epsilon >= 0., "must be >= 0.0");
        check_parameter("TP_Tiny", tp_params.TP_Tiny, tp_params.TP_Tiny >= 0.,
                        "must be >= 0.0");
        check_parameter("TP_Extend_Radius", tp_params.TP_Extend_Radius,
                        tp_params.TP_Extend_Radius >= 0., "must be >= 0.0");
    }
    double tp_offset_plus, tp_offset_minus;
    TP::Parameters tp_params;
#endif
};

#endif /* SIMULATIONPARAMETERS_HPP */
