#ifndef KLEINGORDONLEVEL_HPP_
#define KLEINGORDONLEVEL_HPP_

#include "DefaultLevelFactory.hpp"
#include "DerivedVariables.hpp"
#include "GRAMRLevel.hpp"
#include "InitialConditions.hpp"
#include "KleinGordonRHS.hpp"
#include "Potential.hpp"
#include "VarsTools.hpp"

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
    void specificUpdateODE(amrex::MultiFab &a_soln) override{};

    // to do post each time step on every level
    void specificPostTimeStep() override{};

    //! Error estimation for regridding.
    void errorEst(amrex::TagBoxArray &tags, int clearval, int tagval,
                  amrex::Real time, int n_error_buf = 0,
                  int ngrow = 0) override;

  private:

    KleinGordonLevel &getLevel(int lev)
    {
        return dynamic_cast<KleinGordonLevel &>(parent->getLevel(lev));
    }
};

#endif /* KLEINGORDONLEVEL_HPP_ */
