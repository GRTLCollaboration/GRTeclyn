#ifndef SUPPORTEDWORMHOLELEVEL_HPP_
#define SUPPORTEDWORMHOLELEVEL_HPP_

#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"

class SupportedWormholeLevel : public GRAMRLevel
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
};

#endif /* SUPPORTEDWORMHOLELEVEL_HPP_ */