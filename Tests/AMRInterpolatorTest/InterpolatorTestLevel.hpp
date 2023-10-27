/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INTERPOLATORTESTLEVEL_HPP_
#define INTERPOLATORTESTLEVEL_HPP_

#include "GRAMRLevel.hpp"
// #include "Polynomial.hpp"

class InterpolatorTestLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<InterpolatorTestLevel>;
    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    // initialize data
    void initData() override
    {
        // BoxLoops::loop(Polynomial(m_p.center, m_dx), m_state_new,
        // m_state_new,
        //                FILL_GHOST_CELLS);
    }

    void specificEvalRHS(amrex::MultiFab &a_soln, amrex::MultiFab &a_rhs,
                         const double a_time) override
    {
    }

    void errorEst(amrex::TagBoxArray &tag_box_array, int clearval, int tagval,
                  amrex::Real time, int n_error_buf = 0, int ngrow = 0)
    {
    }
};

#endif /* INTERPOLATORTESTLEVEL_HPP_ */
