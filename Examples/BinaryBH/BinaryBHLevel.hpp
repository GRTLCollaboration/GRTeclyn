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
    friend class DefaultLevelFactory<BinaryBHLevel>;
    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

#if 0
    // xxxxx
    BHAMR &m_bh_amr = dynamic_cast<BHAMR &>(m_gr_amr);
#ifdef USE_TWOPUNCTURES
    TPAMR &m_tp_amr = dynamic_cast<TPAMR &>(m_gr_amr);
#endif /* USE_TWOPUNCTURES */
#endif

    /// Things to do at every full timestep
    ///(might include several substeps, e.g. in RK4)
    virtual void specificAdvance();//xxxxx override;

    /// Initial data calculation
    virtual void initialData();//xxxxx override;

    /// Calculation of the right hand side for the time stepping
    virtual void specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                 const double a_time);//xxxxx override;

    /// Things to do after dt*rhs has been added to the solution
    virtual void specificUpdateODE(GRLevelData &a_soln,
                                   const GRLevelData &a_rhs,
                                   amrex::Real a_dt);//xxxxx override;

    /// Things to do before tagging cells (i.e. filling ghosts)
    virtual void preTagCells();//xxxxx override;

    /// Identify and tag the cells that need higher resolution
    virtual void
    computeTaggingCriterion(amrex::FArrayBox &tagging_criterion,
                            const amrex::FArrayBox &current_state);//xxxxx override;

    // to do post each time step on every level
    virtual void specificPostTimeStep();//xxxxx override;

#ifdef CH_USE_HDF5
    /// Any actions that should happen just before plot files output
    virtual void prePlotLevel() override;
#endif /* CH_USE_HDF5 */
};

#endif /* BINARYBHLEVEL_HPP_ */
