/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP
#define SIMULATIONPARAMETERS_HPP

// General includes
#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"

// Problem specific includes
#include "SphericalExtractionParameters.hpp"

class SimulationParameters : public SimulationParametersBase
{
  public:
    // NOLINTNEXTLINE(readability-identifier-length)
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        read_params(pp);
    }

    /// Read shared parameters
    // NOLINTNEXTLINE(readability-identifier-length)
    void read_params(GRParmParse &pp)
    {

        // Do we want puncture tracking and constraint norm calculation?
        pp.load("puncture_tracking.enabled", puncture_tracking_enabled, true);
        pp.load("puncture_tracking.level", puncture_tracking_level, max_level);
        pp.load("puncture_tracking.initial_coords",
                puncture_tracking_initial_coords,
                {center[0], center[1] - 1.0, center[2], center[0],
                 center[1] + 1.0, center[2]});

        pp.load("fake_bh1_mass", fake_bh1_mass, 0.5);
        pp.load("fake_bh2_mass", fake_bh2_mass, 0.5);

        pp.load("num_points", num_points, 2);
        pp.load("verbosity", verbosity, false);

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
};

#endif /* SIMULATIONPARAMETERS_HPP */
