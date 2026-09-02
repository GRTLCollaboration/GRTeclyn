/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef TEUKOLSKYWAVELEVEL_HPP_
#define TEUKOLSKYWAVELEVEL_HPP_

#include "DefaultLevelBld.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRAmrLevel.hpp"
#include "SixthOrderDerivatives.hpp"
#include "TeukolskyWaveAmr.hpp"

/// Evolution level for a teukolsky wave.
class TeukolskyWaveLevel : public GRAmrLevel
{
  public:
    using GRAmrLevel::GRAmrLevel;

    /// GRAmrLevel only knows about the base GRAmr; the Weyl4 extraction
    /// interpolator lives on the TeukolskyWaveAmr subclass constructed in
    /// Main_TeukolskyWave.cpp, so we need to cast back down to it.
    TeukolskyWaveAmr *get_tw_amr_ptr();

    static void variableSetUp();

    void specific_advance() override;

    void initData() override;

    void specific_eval_rhs(amrex::MultiFab &a_soln, amrex::MultiFab &a_rhs,
                           amrex::Real a_time) override;

    void specific_update_ode(amrex::MultiFab &a_soln) override;

    void specific_post_timestep() override;

    void tag_cells(amrex::TagBoxArray &a_tag_box_array,
                   amrex::Real a_regrid_threshold) final;
};

#endif /* TEUKOLSKYWAVELEVEL_HPP_ */
