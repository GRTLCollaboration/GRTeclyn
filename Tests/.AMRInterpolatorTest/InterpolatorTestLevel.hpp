/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef INTERPOLATORTESTLEVEL_HPP_
#define INTERPOLATORTESTLEVEL_HPP_

#include "BoxLoops.hpp"
#include "GRAMRLevel.hpp"
#include "Polynomial.hpp"
#include "SetValue.hpp"

class InterpolatorTestLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<InterpolatorTestLevel>;
    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    // initialize data
    virtual void initialData()
    {
        BoxLoops::loop(Polynomial(m_p.center, m_dx), m_state_new, m_state_new,
                       FILL_GHOST_CELLS);
    }

    virtual void specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                 const double a_time)
    {
    }

    virtual void computeTaggingCriterion(FArrayBox &tagging_criterion,
                                         const FArrayBox &current_state){};
};

#endif /* INTERPOLATORTESTLEVEL_HPP_ */
