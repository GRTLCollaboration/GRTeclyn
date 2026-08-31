/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef BHAMR_HPP_
#define BHAMR_HPP_

#include "GRAmr.hpp"
#include "ParticleInterpolator.hpp"
#include "PunctureTracker.hpp"

/// A child of Chombo's AMR class to interface with tools which require
/// access to the whole AMR hierarchy, and those of GRAmr
/**
 * This object inherits from GRAmr and adds tools required for BH spacetimes
 */

template <int num_punctures> class BHAmr : public GRAmr
{
  private:
    PunctureTracker<num_punctures> m_puncture_tracker;

  public:

    // example for interpolator object for Psi4 extraction
    static constexpr int weyl_num_components = 2;
    ParticleInterpolator<weyl_num_components>
        m_weyl_interpolator; // interpolator object

    BHAmr(amrex::LevelBld *a_levelbld) : GRAmr(a_levelbld)
    {
        m_puncture_tracker.configure();
        if (m_puncture_tracker.is_enabled())
        {
            m_puncture_tracker.initialize(this);
        }
    }

    void init(amrex::Real a_strt_time, amrex::Real a_stop_time) override
    {
        GRAmr::init(a_strt_time, a_stop_time);

        m_weyl_interpolator.setup(this);
    }

    PunctureTracker<num_punctures> &get_puncture_tracker()
    {
        return m_puncture_tracker;
    }

    [[nodiscard]] bool puncture_tracking_enabled() const
    {
        return m_puncture_tracker.is_enabled();
    }
};

#endif /* BHAMR_HPP_ */
