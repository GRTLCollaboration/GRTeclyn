/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef BINARYBHLEVEL_HPP_
#define BINARYBHLEVEL_HPP_

#include "BHAMR.hpp"
#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"

class BinaryBHLevel : public GRAMRLevel
{
  public:
    static void variableSetUp();

    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    static constexpr int num_punctures = 2;

    BHAMR<num_punctures> *get_bhamr_ptr();

    /// Get a reference to the PunctureTracker object stored by BHAMR
    PunctureTracker<num_punctures> &get_puncture_tracker();

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

    /// Things to do before tagging cells for regridding
    void pre_tag_cells() final;

    /// Tag cells for regridding
    void tag_cells(amrex::TagBoxArray &a_tag_box_array,
                   amrex::Real a_regrid_threshold) final;

    //! Things to do after a restart
    void specific_post_restart() override;

    //! Things to do after init
    void specific_post_init() override;

    //! Things to do after writing a plotfile
    void specific_post_plotfile(const std::string &a_dir,
                                std::ostream & /*a_os*/) override;

    //! Things to do after writing a checkpoint
    void specific_post_checkpoint(const std::string &a_dir,
                                  std::ostream & /*a_os*/) override;
};

#endif /* BINARYBHLEVEL_HPP_ */
