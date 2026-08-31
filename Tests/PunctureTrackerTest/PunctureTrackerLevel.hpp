/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef PUNCTURETRACKERLEVEL_HPP_
#define PUNCTURETRACKERLEVEL_HPP_

#include "BHAmr.hpp"
#include "DefaultLevelBld.hpp"
#include "GRAmrLevel.hpp"

class PunctureTrackerLevel : public GRAmrLevel
{
  public:
    static void variableSetUp();

    // Inherit the contructors from GRAmrLevel
    using GRAmrLevel::GRAmrLevel;

    static constexpr int num_punctures = 2;
    static constexpr std::size_t num_puncture_coords =
        static_cast<std::size_t>(AMREX_SPACEDIM * num_punctures);
    static constexpr amrex::Real shift_y_val = -1.0;

    BHAmr<num_punctures> *get_bh_amr_ptr();

    /// Get a reference to the PunctureTracker object stored by BHAmr
    PunctureTracker<num_punctures> &get_puncture_tracker();

    /// Initial data calculation
    void initData() override;

    /// Calculation of the right hand side for the time stepping
    void specific_eval_rhs(amrex::MultiFab &a_soln, amrex::MultiFab &a_rhs,
                           const amrex::Real a_time) override;

    // to do post each time step on every level
    void specific_post_timestep() override;

    /// Tag cells for regridding
    void tag_cells(amrex::TagBoxArray &a_tag_box_array,
                   amrex::Real a_regrid_threshold) final;

    //! Things to do after a restart
    void specific_post_restart() override;

    //! Things to do after init
    void specific_post_init() override;

    //! Things to do after regridding
    void specific_post_regrid(int a_lbase, int a_new_finest) override;

    //! Things to do after writing a checkpoint
    void specific_post_checkpoint(const std::string &a_chk_dir,
                                  std::ostream & /*a_os*/) override;

  private:
    void check_puncture_tagging();
};

#endif /* PUNCTURETRACKERLEVEL_HPP_ */
