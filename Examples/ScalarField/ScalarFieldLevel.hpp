/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SCALARFIELDLEVEL_HPP_
#define SCALARFIELDLEVEL_HPP_

#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"
// Problem specific includes
#include "Potential.hpp"
#include "ScalarField.hpp"

//!  A class for the evolution of a scalar field, minimally coupled to gravity
/*!
     The class takes some initial data for a scalar field (variables phi and Pi)
     and evolves it using the CCZ4 equations. It is possible to specify an
   initial period of relaxation for the conformal factor chi, for non analytic
   initial conditions (for example, a general field configuration at a moment of
   time symmetry assuming conformal flatness). \sa MatterCCZ4(),
   ConstraintsMatter(), ScalarField(), RelaxationChi()
*/
class ScalarFieldLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<ScalarFieldLevel>;

  public:

    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    static void variableSetUp();

    // Typedef for scalar field
    typedef ScalarField<Potential> ScalarFieldWithPotential;

    //! Things to do at the end of the advance step, after RK4 calculation
    void specificAdvance() override;

    //! Initialize data for the field and metric variables

    void initData() override;

    //! RHS routines used at each RK4 step
    void specificEvalRHS(amrex::MultiFab &a_soln, amrex::MultiFab &a_rhs,
                         const double a_time) override;

    //! Things to do in UpdateODE step, after soln + rhs update
    void specificUpdateODE(amrex::MultiFab &a_soln) override;

    /// Things to do before tagging cells (i.e. filling ghosts)
    void preTagCells();

    //! Tell GRTeclyn how to tag cells for regridding
    void errorEst(amrex::TagBoxArray &tag_box_array, int clearval, int tagval,
                  amrex::Real time, int n_error_buf = 0, int ngrow = 0) final;
};

#endif /* SCALARFIELDLEVEL_HPP_ */
