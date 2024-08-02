/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
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

    BHAMR *get_bhamr_ptr();

    /// Things to do at every full timestep
    ///(might include several substeps, e.g. in RK4)
    void specificAdvance() override;

    /// Initial data calculation
    void initData() override;

    /// Calculation of the right hand side for the time stepping
    void specificEvalRHS(amrex::MultiFab &a_soln, amrex::MultiFab &a_rhs,
                         const double a_time) override;

    /// Things to do after dt*rhs has been added to the solution
    void specificUpdateODE(amrex::MultiFab &a_soln) override;

    // to do post each time step on every level
    void specificPostTimeStep() override;

    void errorEst(amrex::TagBoxArray &tag_box_array, int clearval, int tagval,
                  amrex::Real time, int n_error_buf = 0, int ngrow = 0) final;

    //! Things to do after a restart
    void specific_post_restart() override;

    //! Things to do after init
    void specific_post_init() override;

    //! Things to do after writing a checkpoint
    void specificPostCheckpoint(const std::string &a_dir,
                                std::ostream & /*a_os*/) override;

#ifdef AMREX_USE_HDF5
    /// Any actions that should happen just before plot files output
    virtual void prePlotLevel() override;
#endif /* AMREX_USE_HDF5 */

  private:
    void restart_punctures();
};

#endif /* BINARYBHLEVEL_HPP_ */
