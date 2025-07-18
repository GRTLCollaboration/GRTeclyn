#ifndef KLEINGORDONLEVEL_HPP_
#define KLEINGORDONLEVEL_HPP_

#include "DefaultLevelFactory.hpp"
#include "DerivedVariables.hpp"
#include "GRAMRLevel.hpp"
#include "KleinGordonRHS.hpp"
#include "VarsTools.hpp"

#include <variant>

class KleinGordonLevel : public GRAMRLevel
{
  public:
    using GRAMRLevel::GRAMRLevel;

    //! Define data descriptors.
    static void variableSetUp();
    //    static void variableCleanUp ();

    //! Initialize data at problem start-up.
    void initData() override;

    //! Advance this level for one step

    void specificEvalRHS(amrex::MultiFab &a_soln, amrex::MultiFab &a_rhs,
                         const double a_time) override;

    void specificAdvance() override;

    /// Things to do after dt*rhs has been added to the solution
    void specificUpdateODE(amrex::MultiFab &a_soln) override {};

    // to do post each time step on every level
    void specificPostTimeStep() override {};

    //! Error estimation for regridding.
    void tag_cells(amrex::TagBoxArray &tags,
                   amrex::Real a_regrid_threshold) override;

  private:

    KleinGordonLevel &getLevel(int lev)
    {
        return dynamic_cast<KleinGordonLevel &>(parent->getLevel(lev));
    }
};

#endif /* KLEINGORDONLEVEL_HPP_ */
