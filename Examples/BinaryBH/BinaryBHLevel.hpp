/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BINARYBHLEVEL_HPP_
#define BINARYBHLEVEL_HPP_

#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"
// TPAMR.hpp includes BHAMR.hpp
#include "TPAMR.hpp"

class BinaryBHLevel : public GRAMRLevel
{
  public:

    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    /// Things to do at every full timestep
    ///(might include several substeps, e.g. in RK4)
    void specificAdvance() override;

    /// Initial data calculation
    virtual void initData() override;

    /// Calculation of the right hand side for the time stepping
    virtual void specificEvalRHS(amrex::MultiFab &a_soln,
                                 amrex::MultiFab &a_rhs,
                                 const double a_time) override;

    /// Things to do after dt*rhs has been added to the solution
    virtual void specificUpdateODE(amrex::MultiFab &a_soln) override;

    // to do post each time step on every level
    virtual void specificPostTimeStep() override;

    virtual void errorEst(amrex::TagBoxArray &tb, int clearval, int tagval,
                          amrex::Real time, int n_error_buf = 0,
                          int ngrow = 0) override final;

#ifdef AMREX_USE_HDF5
    /// Any actions that should happen just before plot files output
    virtual void prePlotLevel() override;
#endif /* AMREX_USE_HDF5 */
};

#endif /* BINARYBHLEVEL_HPP_ */
