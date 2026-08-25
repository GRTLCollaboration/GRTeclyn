/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SCALARFIELDAMR_HPP_
#define SCALARFIELDAMR_HPP_

#include "FilesystemTools.hpp"
#include "GRAMR.hpp"
#include "GRParmParse.hpp"
#include "ParticleInterpolator.hpp"

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

        GRParmParse line_pp("line_extraction");
        bool line_extraction_enabled{};
        line_pp.get("enabled", line_extraction_enabled);
        if (line_extraction_enabled)
        {
            GRParmParse grteclyn_pp("grteclyn");
            std::string output_path{};
            grteclyn_pp.get("output_path", output_path);
            std::string output_subpath{};
            line_pp.get("output_subpath", output_subpath);
            FilesystemTools::ensure_directory_exists(output_path + "/" +
                                                     output_subpath);

            phi_interpolator.setup(this);
            rho_interpolator.setup(this);
        }
    }
};

#endif /* SCALARFIELDAMR_HPP_ */
