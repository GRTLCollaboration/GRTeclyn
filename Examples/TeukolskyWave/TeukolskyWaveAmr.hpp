/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef TEUKOLSKYWAVEAMR_HPP_
#define TEUKOLSKYWAVEAMR_HPP_

#include "GRAmr.hpp"
#include "ParticleInterpolator.hpp"

/// A child of GRAmr that adds the particle interpolator needed for the Weyl4
/// extraction used by TeukolskyWaveLevel::specific_post_timestep().
class TeukolskyWaveAmr : public GRAmr
{
  public:
    static constexpr int weyl_num_components = 2;
    ParticleInterpolator<weyl_num_components>
        m_weyl_interpolator; // interpolator object

    using GRAmr::GRAmr;

    void init(amrex::Real a_strt_time, amrex::Real a_stop_time) override
    {
        GRAmr::init(a_strt_time, a_stop_time);
        m_weyl_interpolator.setup(this);
    }
};

#endif /* TEUKOLSKYWAVEAMR_HPP_ */
