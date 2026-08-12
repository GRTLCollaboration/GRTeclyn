/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SCALARFIELDLEVEL_HPP_
#define SCALARFIELDLEVEL_HPP_

#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"
#include "Potential.hpp"
#include "ScalarField.hpp"

/// Evolution level for a real scalar field minimally coupled to gravity.
class ScalarFieldLevel : public GRAMRLevel
{
  public:
    using GRAMRLevel::GRAMRLevel;

    using ScalarFieldWithPotential = ScalarField<Potential>;

    static void variableSetUp();

    void specificAdvance() override;

    void initData() override;

    void specificEvalRHS(amrex::MultiFab &a_soln, amrex::MultiFab &a_rhs,
                         double a_time) override;

    void specificUpdateODE(amrex::MultiFab &a_soln) override;

    void tag_cells(amrex::TagBoxArray &a_tag_box_array,
                   amrex::Real a_regrid_threshold) final;
};

#endif /* SCALARFIELDLEVEL_HPP_ */
