/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef BHAMR_HPP_
#define BHAMR_HPP_

#include "GRAMR.hpp"
#include "PunctureTracker.hpp"

/// A child of Chombo's AMR class to interface with tools which require
/// access to the whole AMR hierarchy, and those of GRAMR
/**
 * This object inherits from GRAMR and adds tools required for BH spacetimes
 */
class BHAMR : public GRAMR
{
  public:
    PunctureTracker m_puncture_tracker;

    BHAMR(amrex::LevelBld *a_levelbld) : GRAMR(a_levelbld) {}
};

#endif /* BHAMR_HPP_ */
