#ifndef KLEINGORDONLEVEL_HPP_
#define KLEINGORDONLEVEL_HPP_

#include "DefaultLevelFactory.hpp"
#include "DerivedVariables.hpp"
#include "GRAMRLevel.hpp"
#include "InitialConditions.H"
#include "Potential.H"

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
    // could also use GRAMRLevel::advance
    // but would have to also define specificEvalRHS and specificAdvance and
    // specificUpdateODE
    //  amrex::Real advance(amrex::Real time, amrex::Real dt, int iteration,
    //                      int ncycle) override;

    void specificEvalRHS(amrex::MultiFab &a_soln, amrex::MultiFab &a_rhs,
                         const double a_time) override;

    void specificAdvance() override;

    /// Things to do after dt*rhs has been added to the solution
    void specificUpdateODE(amrex::MultiFab &a_soln) override{};

    // to do post each time step on every level
    void specificPostTimeStep() override{};

    //! Error estimation for regridding.
    void errorEst(amrex::TagBoxArray &tb, int clearval, int tagval,
                  amrex::Real time, int n_error_buf = 0,
                  int ngrow = 0) override;

  private:

    static amrex::Vector<std::string> diagnostics; // this is for error checking

    //! Get AmrLevelWave
    KleinGordonLevel &getLevel(int lev)
    {
        return static_cast<KleinGordonLevel &>(parent->getLevel(lev));
    }
};

#endif /* KLEINGORDONLEVEL_HPP_ */
