/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP
#define SIMULATIONPARAMETERS_HPP

// General includes
#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"

// Problem specific includes:
#include "ArrayTools.hpp"
#include "BoostedBHInitialData.hpp"
#ifdef USE_TWOPUNCTURES
#include "TwoPuncturesInitialData.hpp"
#endif

class SimulationParameters : public SimulationParametersBase
{
  public:
    // NOLINTNEXTLINE(readability-identifier-length)
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        read_shared_params(pp);

        read_bh_params(pp);

        check_params();
    }

    /// Read shared parameters
    // NOLINTNEXTLINE(readability-identifier-length)
    void read_shared_params(GRParmParse &pp)
    {
        // Do we want puncture tracking and constraint norm calculation?
        pp.load("puncture_tracking.enabled", puncture_tracking_enabled, false);
        pp.load("puncture_tracking.level", puncture_tracking_level, max_level);
        pp.load("puncture_tracking.writeout_level",
                puncture_tracking_writeout_level, 0);
        pp.load("calculate_constraint_norms", calculate_constraint_norms,
                false);
    }

    /// Read BH parameters
    // NOLINTNEXTLINE(readability-*)
    void read_bh_params(GRParmParse &pp)
    {
#ifndef USE_TWOPUNCTURES
        // Initial data
        pp.load("massA", bh1_params.mass);
        pp.load("momentumA", bh1_params.momentum);
        pp.load("massB", bh2_params.mass);
        pp.load("momentumB", bh2_params.momentum);

        // Get the centers of the BHs either explicitly or as
        // an offset (not both, or they will be offset from center
        // provided)
        std::array<double, AMREX_SPACEDIM> centerA{};
        std::array<double, AMREX_SPACEDIM> centerB{};
        std::array<double, AMREX_SPACEDIM> offsetA{};
        std::array<double, AMREX_SPACEDIM> offsetB{};
        pp.load("centerA", centerA, center);
        pp.load("centerB", centerB, center);
        pp.load("offsetA", offsetA, {0.0, 0.0, 0.0});
        pp.load("offsetB", offsetB, {0.0, 0.0, 0.0});
        FOR (idir)
        {
            bh1_params.center[idir] = centerA[idir] + offsetA[idir];
            bh2_params.center[idir] = centerB[idir] + offsetB[idir];
        }
#else
        // Workaround for setting normal BH parameters when using TwoPunctures
        // as these are used for e.g. puncture tracking and tagging.
        // TODO: Make this less hacky when refactoring parameters
        TwoPuncturesInitialData::read_parameters();
        const auto &two_punctures =
            TwoPuncturesInitialData::get_two_punctures();

        if (two_punctures.give_bare_mass)
        {
            bh1_params.mass = two_punctures.par_m_minus;
            bh2_params.mass = two_punctures.par_m_plus;
        }
        else
        {
            bh1_params.mass = two_punctures.target_M_minus;
            bh2_params.mass = two_punctures.target_M_plus;
        }

        FOR (idir)
        {
            bh1_params.momentum[idir] = two_punctures.par_P_minus[idir];
            bh2_params.momentum[idir] = two_punctures.par_P_plus[idir];
        }

        bh1_params.center = center;
        bh2_params.center = center;
        int offset_dir    = (two_punctures.swap_xz) ? 2 : 0;
        bh1_params.center[offset_dir] +=
            two_punctures.center_offset[offset_dir] - two_punctures.par_b;
        bh2_params.center[offset_dir] +=
            two_punctures.center_offset[offset_dir] + two_punctures.par_b;
#endif
    }

    void check_params()
    {
#if 0
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
#endif

        warn_parameter("massA", bh1_params.mass, bh1_params.mass >= 0,
                       "should be >= 0");
        warn_parameter("massB", bh2_params.mass, bh2_params.mass >= 0,
                       "should be >= 0");
        warn_array_parameter(
            "momentumA", bh1_params.momentum,
            std::sqrt(ArrayTools::norm2(bh1_params.momentum)) <
                0.3 * bh1_params.mass,
            "approximation used for boosted BH only valid for small boosts");
        warn_array_parameter(
            "momentumB", bh2_params.momentum,
            std::sqrt(ArrayTools::norm2(bh2_params.momentum)) <
                0.3 * bh1_params.mass,
            "approximation used for boosted BH only valid for small boosts");
        FOR (idir)
        {
            std::string nameA = "centerA[" + std::to_string(idir) + "]";
            std::string nameB = "centerB[" + std::to_string(idir) + "]";
            // NOLINTBEGIN(cppcoreguidelines-init-variables)
            double center_A_dir = bh1_params.center[idir];
            double center_B_dir = bh2_params.center[idir];
            // NOLINTEND(cppcoreguidelines-init-variables)
            warn_parameter(nameA, center_A_dir,
                           (center_A_dir >= 0.0) &&
                               (center_A_dir <= (ivN[idir] + 1) * coarsest_dx),
                           "should be within the computational domain");
            warn_parameter(nameB, center_B_dir,
                           (center_B_dir >= 0.0) &&
                               (center_B_dir <= (ivN[idir] + 1) * coarsest_dx),
                           "should be within the computational domain");
        }
        check_parameter("puncture_tracking_level", puncture_tracking_level,
                        (puncture_tracking_level >= 0) &&
                            (puncture_tracking_level <= max_level),
                        "must be between 0 and max_level (inclusive)");
    }

    bool puncture_tracking_enabled{};
    int puncture_tracking_level{};
    int puncture_tracking_writeout_level{};
    bool calculate_constraint_norms{};

    // Collection of parameters necessary for initial conditions
    // Set these even in the case of TwoPunctures as they are used elsewhere
    // e.g. for puncture tracking/tagging
    BoostedBHInitialData::params_t bh2_params{};
    BoostedBHInitialData::params_t bh1_params{};
};

#endif /* SIMULATIONPARAMETERS_HPP */
