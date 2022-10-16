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
    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

#if 0
    friend class DefaultLevelFactory<BinaryBHLevel>;
    // xxxxx
    BHAMR &m_bh_amr = dynamic_cast<BHAMR &>(m_gr_amr);
#ifdef USE_TWOPUNCTURES
    TPAMR &m_tp_amr = dynamic_cast<TPAMR &>(m_gr_amr);
#endif /* USE_TWOPUNCTURES */
#endif

    /// Things to do at every full timestep
    ///(might include several substeps, e.g. in RK4)
    virtual void specificAdvance() override;

    /// Initial data calculation
    virtual void initData() override;

    /// Calculation of the right hand side for the time stepping
    virtual void specificEvalRHS(amrex::MultiFab& a_soln,
                                 amrex::MultiFab& a_rhs,
                                 const double a_time) override;

    /// Things to do after dt*rhs has been added to the solution
    virtual void specificUpdateODE(amrex::MultiFab& a_soln) override;

    /// Things to do before tagging cells (i.e. filling ghosts)
    virtual void preTagCells();//xxxxx override;

    /// Identify and tag the cells that need higher resolution
    virtual void
    computeTaggingCriterion(amrex::FArrayBox &tagging_criterion,
                            const amrex::FArrayBox &current_state);//xxxxx override;

    // to do post each time step on every level
    virtual void specificPostTimeStep() override;

#ifdef AMREX_USE_HDF5
    /// Any actions that should happen just before plot files output
    virtual void prePlotLevel() override;
#endif /* AMREX_USE_HDF5 */
};

#endif /* BINARYBHLEVEL_HPP_ */
