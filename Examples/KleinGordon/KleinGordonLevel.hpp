#ifndef KLEINGORDONLEVEL_HPP_
#define KLEINGORDONLEVEL_HPP_

#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"
#include "InitialConditions.H"
#include "KleinGordonDerive.hpp"
#include "Potential.H"

class KleinGordonLevel : public GRAMRLevel
{
  public:
    using GRAMRLevel::GRAMRLevel;

    //! Define data descriptors.
    //    static void variableSetUp ();
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

    void derive(const std::string &name, amrex::Real /*time*/,
                amrex::MultiFab &multifab, int dcomp) override;

    //    Need to be public for CUDA
    //     void computeRHS (amrex::MultiFab& dSdt, amrex::MultiFab const& S);

  private:

    //   static constexpr int nghost = 2; // Two ghost cells needed

    //   static int verbose;
    //   static int rk_order;
    //   static amrex::Real cfl;
    // //    static double ampl;
    // //    static double width;
    //   static amrex::Vector<float> ampl;
    //   static amrex::Vector<float> width;
    //   static int nfields;
    //   static amrex::Real scalar_mass;
    //   static amrex::Real k_r;
    //   static amrex::Real alpha;

    //   static int ncomp;
    static amrex::Vector<std::string>
        plot_constraints; // this is for error checking
    // //  static constexpr int ncomp;//  = nfields*2; // Two variables: u & v

    //! Get AmrLevelWave
    KleinGordonLevel &getLevel(int lev)
    {
        return static_cast<KleinGordonLevel &>(parent->getLevel(lev));
    }
};

#endif /* KLEINGORDONLEVEL_HPP_ */
