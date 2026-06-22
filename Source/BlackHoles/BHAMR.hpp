/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef BHAMR_HPP_
#define BHAMR_HPP_

#include "GRAMR.hpp"
#include "ParticleInterpolator.hpp"
#include "PunctureTracker.hpp"

#include <AMReX_ParmParse.H>

/// A child of Chombo's AMR class to interface with tools which require
/// access to the whole AMR hierarchy, and those of GRAMR
/**
 * This object inherits from GRAMR and adds tools required for BH spacetimes
 */

template <int num_punctures> class BHAMR : public GRAMR
{
  private:
    PunctureTracker<num_punctures> m_puncture_tracker;

  public:

    // example for interpolator object for Psi4 extraction
    static constexpr int weyl_num_components = 2;
    ParticleInterpolator<weyl_num_components>
        m_weyl_interpolator; // interpolator object

    BHAMR(amrex::LevelBld *a_levelbld) : GRAMR(a_levelbld)
    {
        amrex::ParmParse puncture_tracking_pp("puncture_tracking");
        bool puncture_tracking_enabled = false; // default

        puncture_tracking_pp.query("enabled", puncture_tracking_enabled);
        if (puncture_tracking_enabled)
        {
            m_puncture_tracker.initialize(this);
        }
    }

    void init(amrex::Real a_strt_time, amrex::Real a_stop_time) override
    {
        GRAMR::init(a_strt_time, a_stop_time);

        const auto &params = get_simulation_parameters();
        m_weyl_interpolator.setup(this, params.boundary_params, true);
    }

    PunctureTracker<num_punctures> &get_puncture_tracker()
    {
        return m_puncture_tracker;
    }
};

#endif /* BHAMR_HPP_ */
