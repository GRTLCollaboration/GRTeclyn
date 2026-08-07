/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "BaseParameterChecker.hpp"
#include "GRParmParse.hpp"

// Problem specific includes
#include "SphericalExtraction.hpp"

class SimulationParameters : public BaseParameterChecker
{
  public:
    SimulationParameters() : BaseParameterChecker() { check_params(pp); }

  private:
    void check_params()
    {
        GRParmParse pp;
        // Extraction params
        pp.load("num_extraction_radii",
                extraction_params_lo.num_extraction_radii(), 1);

        if (pp.contains("extraction_radii"))
        {
            pp.load("extraction_radii", extraction_params_lo.extraction_radii(),
                    extraction_params_lo.num_extraction_radii());
        }
        else
        {
            pp.load("extraction_radius",
                    extraction_params_lo.extraction_radii(), 1, L / 4);
        }
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

  public:
    SphericalExtraction::params_t extraction_params_lo;
    int es, el, em; // spherical harmonic params
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
