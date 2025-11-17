/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef BHAMR_HPP_
#define BHAMR_HPP_

#include "GRAMR.hpp"
#include "InterpolationQueryParticle.hpp"
#include "ParticleInterpolators.hpp"
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
    InterpolationQueryParticle *m_query =
        nullptr; // query used for interpolation
    bool m_query_populated =
        false; // flag to identify whether the query has been populated: for
               // particles that will be fixed we want to do this only once
               // here!

  public:

    ParticleInterpolators<1> *m_weyl_interpolator =
        nullptr; // weyl interpolator (used as chi interpolator in this example)

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

    PunctureTracker<num_punctures> &get_puncture_tracker()
    {
        return m_puncture_tracker;
    }

    // set weyl interpolator
    void set_interpolator(ParticleInterpolators<1> *a_interpolator)
    {
        AMREX_ASSERT(a_interpolator != nullptr);

        m_weyl_interpolator = a_interpolator;
        m_weyl_interpolator->set_gramr_ptr(this);
    }

    // set query
    void set_query(InterpolationQueryParticle &q) override
    {
        if (m_query != &q)
        { // only if the query object changed
            amrex::Print() << "Query object changed \n";
            m_query           = &q;
            m_query_populated = false;
        }
    }

    // populate only once
    void ensure_query_populated() override
    {
        AMREX_ASSERT(m_query);
        if (!m_query_populated)
        {
            amrex::Print() << "Query being populated \n";
            m_weyl_interpolator->populate_from_query(*m_query);
            m_query_populated = true;
        }
    }

    // access to a cached query
    InterpolationQueryParticle *query() override
    {
        return m_query ? &*m_query : nullptr;
    }
};

#endif /* BHAMR_HPP_ */
