/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef INTERPOLATORTESTLEVEL_HPP_
#define INTERPOLATORTESTLEVEL_HPP_

#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"
#include "PolynomialTest.hpp"

// We basically need this to have a valid AMR hierarchy

class InterpolatorTestLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<InterpolatorTestLevel>;

    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    static void variableSetUp()
    {

        stateVariableSetUp();
        PolynomialTest::set_up(State_Type);
    }

    // initialize data
    void initData()
    {
        // We do not need initial data I think. I precompute the polynomial
    }

    void specificEvalRHS(amrex::MultiFab &a_soln, amrex::MultiFab &a_rhs,
                         const double a_time)
    {
    }

    void tag_cells(amrex::TagBoxArray &a_tag_box_array,
                   amrex::Real a_regrid_threshold) final
    {
    }
};

#endif /* INTERPOLATORTESTLEVEL_HPP_ */