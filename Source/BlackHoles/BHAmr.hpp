/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef BHAMR_HPP_
#define BHAMR_HPP_

#include "AHFinder.hpp"
#include "GRAmr.hpp"
#include "ParticleInterpolator.hpp"
#include "PunctureTracker.hpp"

#include <memory>
#include <utility>
#include <vector>

/// A child of Chombo's AMR class to interface with tools which require
/// access to the whole AMR hierarchy, and those of GRAmr
/**
 * This object inherits from GRAmr and adds tools required for BH spacetimes
 */

template <int num_punctures> class BHAmr : public GRAmr
{
  private:
    PunctureTracker<num_punctures> m_puncture_tracker;

    // One apparent horizon finder per puncture
    std::vector<std::unique_ptr<AHFinder>> m_ah_finders;

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

    [[nodiscard]] bool ah_finder_enabled() const
    {
        return !m_ah_finders.empty();
    }

    // Take ownership of an apparent horizon finder and initialise it
    AHFinder &add_ah_finder(std::unique_ptr<AHFinder> a_ah_finder)
    {
        m_ah_finders.push_back(std::move(a_ah_finder));
        m_ah_finders.back()->init(this);
        return *m_ah_finders.back();
    }

    [[nodiscard]] int num_ah_finders() const
    {
        return static_cast<int>(m_ah_finders.size());
    }

    // The i-th finder (0-based), if ah_finder_enabled()
    AHFinder &get_ah_finder(int i) { return *m_ah_finders.at(i); }
};

#endif /* BHAMR_HPP_ */
