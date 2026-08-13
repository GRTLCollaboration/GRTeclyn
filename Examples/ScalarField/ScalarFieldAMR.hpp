/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SCALARFIELDAMR_HPP_
#define SCALARFIELDAMR_HPP_

#include "FilesystemTools.hpp"
#include "GRAMR.hpp"
#include "ParticleInterpolator.hpp"
#include "SimulationParameters.hpp"

//! AMR hierarchy carrying the interpolators used for line extraction.
class ScalarFieldAMR : public GRAMR
{
  public:
    ParticleInterpolator<1> phi_interpolator;
    ParticleInterpolator<1> rho_interpolator;

    explicit ScalarFieldAMR(amrex::LevelBld *a_level_bld) : GRAMR(a_level_bld)
    {
    }

    void init(amrex::Real a_start_time, amrex::Real a_stop_time) override
    {
        GRAMR::init(a_start_time, a_stop_time);

        const auto &params = get_simulation_parameters();
        if (params.activate_line_extraction)
        {
            if (!params.data_path.empty() &&
                !FilesystemTools::directory_exists(params.data_path))
            {
                FilesystemTools::mkdir_recursive(params.data_path);
            }

            phi_interpolator.setup(this, params.boundary_params);
            rho_interpolator.setup(this, params.boundary_params);
        }
    }
};

#endif /* SCALARFIELDAMR_HPP_ */
