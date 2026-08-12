/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP
#define SIMULATIONPARAMETERS_HPP

// General includes
#include "BaseParameterChecker.hpp"
#include "GRParmParse.hpp"
#include "PunctureTracker.hpp"
#include "CCZ4RHS.hpp"

class SimulationParameters
{
  public:
    // NOLINTNEXTLINE(readability-identifier-length)
    SimulationParameters() = delete;

    /// Read shared parameters
    // NOLINTNEXTLINE(readability-identifier-length)
    static void check_params()
    {
        BaseParameterChecker::check_params();

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

        bool puncture_tracking_enabled{false};
        pp.queryAdd("puncture_tracking.enabled", puncture_tracking_enabled);
        if (puncture_tracking_enabled)
        {
            puncture_tracker_params_t::check_params();
        }

        amrex::Real fake_bh1_mass = 0.5;
        amrex::Real fake_bh2_mass = 0.5;

        pp.queryAdd("fake_bh1_mass", fake_bh1_mass);
        pp.queryAdd("fake_bh2_mass", fake_bh2_mass);

        int num_points = 2;
        pp.queryAdd("num_points", num_points);

        bool verbosity = true;
        pp.queryAdd("verbosity", verbosity);
    }
/*
        // Extraction params
        pp.load("num_extraction_radii",
                extraction_params_lo.num_extraction_radii(), 1);

        std::vector<double> extraction_radii_stdvect;
        if (pp.contains("extraction_radii"))
        {
            pp.load("extraction_radii", extraction_radii_stdvect,
                    extraction_params_lo.num_extraction_radii());
        }
        else
        {
            pp.load("extraction_radius", extraction_radii_stdvect, 1, 0.1);
        }
        extraction_params_lo.extraction_radii().resize(
            extraction_params_lo.num_extraction_radii());
        std::copy(extraction_radii_stdvect.begin(),
                  extraction_radii_stdvect.end(),
                  extraction_params_lo.extraction_radii().begin());

        pp.load("num_points_phi_lo", extraction_params_lo.num_points_phi(), 8);
        pp.load("num_points_theta_lo", extraction_params_lo.num_points_theta(),
                17);
        pp.load("extraction_center", extraction_params_lo.center, center);
        pp.load("write_extraction", extraction_params_lo.write_extraction,
                false);

        pp.load("es", es, 0);
        pp.load("el", el, 2);
        pp.load("em", em, 0);
    }

    bool puncture_tracking_enabled{};
    int puncture_tracking_level{};
    std::array<amrex::Real, AMREX_SPACEDIM * 2UL>
        puncture_tracking_initial_coords{};

    amrex::Real fake_bh1_mass{};
    amrex::Real fake_bh2_mass{};

    // For ParticleInterpolator Test
    int num_points{};
    bool verbosity{};

    // For SphericalExtraction Test
    spherical_extraction_params_t extraction_params_lo;
    int es{}, el{}, em{}; // spherical harmonic params
*/
};

#endif /* SIMULATIONPARAMETERS_HPP */
