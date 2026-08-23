#ifndef RADIALRECIPELEVEL_HPP_
#define RADIALRECIPELEVEL_HPP_

#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"

class RadialRecipeLevel : public GRAMRLevel
{
  public:
    static void variableSetUp();

    using GRAMRLevel::GRAMRLevel;

    void specificAdvance() override;
    void initData() override;
    void specificEvalRHS(amrex::MultiFab &a_soln, amrex::MultiFab &a_rhs,
                         const double a_time) override;
    void specificUpdateODE(amrex::MultiFab &a_soln) override;
    void pre_tag_cells() final;
    void tag_cells(amrex::TagBoxArray &a_tag_box_array,
                   amrex::Real a_regrid_threshold) final;
    void specificPostTimeStep() override;

    //! Persist / restore the recentring box's odometer alongside a checkpoint.
    //! Without it a restarted run resumes with the trajectory silently reset to
    //! zero, which produces a plot that looks plausible and is wrong.
    void specific_pre_checkpoint(const std::string &a_dir,
                                 std::ostream &a_os) override;
    void specific_post_restart() override;

  private:
    //! One-time setup of the recentring box, from the parameters and this
    //! level's grid.  Idempotent: safe to call from several entry points.
    void configure_recentring_box();
};

#endif /* RADIALRECIPELEVEL_HPP_ */
