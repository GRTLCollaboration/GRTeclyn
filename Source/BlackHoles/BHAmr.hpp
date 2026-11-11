/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef BHAMR_HPP_
#define BHAMR_HPP_

#include "AHFinder.hpp"
#include "GRAmr.hpp"
#include "GRParmParse.hpp"
#include "ParticleInterpolator.hpp"
#include "PunctureTracker.hpp"

#include <array>
#include <memory>
#include <string>
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

    void setup_ah_finders()
    {
        GRParmParse ah_finder_pp("ah_finder");

        bool enabled{};
        ah_finder_pp.get("enabled", enabled);
        if (!enabled)
        {
            return;
        }

        int num_particles{};
        ah_finder_pp.get("num_particles", num_particles);

        std::vector<amrex::Real> centers;
        std::vector<amrex::Real> guess_radii;
        ah_finder_pp.getarr("centers", centers);
        ah_finder_pp.getarr("guess_radii", guess_radii);

        if (centers.size() !=
            static_cast<std::size_t>(AMREX_SPACEDIM * num_punctures))
        {
            ah_finder_pp.error("centers", "must have " +
                                              std::to_string(AMREX_SPACEDIM) +
                                              " values per puncture");
        }
        if (guess_radii.size() != static_cast<std::size_t>(num_punctures))
        {
            ah_finder_pp.error("guess_radii", "must have one value per "
                                              "puncture");
        }

        for (int ibh = 0; ibh < num_punctures; ++ibh)
        {
            const std::array<double, AMREX_SPACEDIM> center{
                centers[AMREX_SPACEDIM * ibh + 0],
                centers[AMREX_SPACEDIM * ibh + 1],
                centers[AMREX_SPACEDIM * ibh + 2]};

            m_ah_finders.push_back(std::make_unique<AHFinder>(
                num_particles, center, guess_radii[ibh],
                "ah" + std::to_string(ibh + 1)));
            m_ah_finders.back()->init(this);
        }
    }

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

        setup_ah_finders();
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

    // The finder for puncture i (0-based), if ah_finder_enabled()
    AHFinder &get_ah_finder(int i) { return *m_ah_finders.at(i); }
};

#endif /* BHAMR_HPP_ */
