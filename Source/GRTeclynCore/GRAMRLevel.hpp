/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef GRAMRLEVEL_HPP_
#define GRAMRLEVEL_HPP_

// Other includes
#include "BoundaryConditions.hpp"
#include "GRAMR.hpp"
// xxxxx#include "InterpSource.hpp"
#include "SimulationParameters.hpp"
#include "StateVariables.hpp" // need NUM_VARS

#include <AMReX_AmrLevel.H>

#include <fstream>
#include <limits>
#include <sys/time.h>

// NOLINTNEXTLINE(cppcoreguidelines-special-member-functions)
class GRAMRLevel : public amrex::AmrLevel
{
  public:
    /**
     * \brief Set up the state variables from StateVariables.hpp.
     * This should be called by the child's variableSetUp().
     */
    static void stateVariableSetUp();

    static void variableCleanUp();

    GRAMRLevel();

    GRAMRLevel(amrex::Amr &papa, int lev, const amrex::Geometry &geom,
               const amrex::BoxArray &box_array,
               const amrex::DistributionMapping &distribution_mapping,
               amrex::Real time);

    ~GRAMRLevel() override;

    GRAMR *get_gramr_ptr();

    /**
     * \brief Compute the initial time step.
     */
    void computeInitialDt(int finest_level, int sub_cycle,
                          amrex::Vector<int> &n_cycle,
                          const amrex::Vector<amrex::IntVect> &ref_ratio,
                          amrex::Vector<amrex::Real> &dt_level,
                          amrex::Real stop_time) override;
    /**
     * \brief Compute the next time step.
     */
    void computeNewDt(int finest_level, int sub_cycle,
                      amrex::Vector<int> &n_cycle,
                      const amrex::Vector<amrex::IntVect> &ref_ratio,
                      amrex::Vector<amrex::Real> &dt_min,
                      amrex::Vector<amrex::Real> &dt_level,
                      amrex::Real stop_time, int post_regrid_flag) override;
    /**
     * \brief Do an integration step on this level.  Returns maximum safe
     * time step.  This is a pure virtual function and hence MUST
     * be implemented by derived classes.
     */
    amrex::Real advance(amrex::Real time, amrex::Real dt, int iteration,
                        int ncycle) override;

    /**
     * \brief Contains operations to be done after a timestep.  This is a
     * pure virtual function and hence MUST be implemented by derived
     * classes.
     */
    void post_timestep(int iteration) override;
    /**
     * \brief Operations to be done after regridding
     * This is a pure virtual function and hence MUST be
     * implemented by derived classes.
     */
    void post_regrid(int lbase, int new_finest) override;
    /**
     * \brief Operations to be done after initialization.
     * This is a pure virtual function and hence MUST be
     * implemented by derived classes.
     */
    void post_init(amrex::Real stop_time) override;
    /**
     * \brief Operations to be done after restart.
     */
    void post_restart() override;
    /**
     * \brief Init data on this level from another AmrLevel (during regrid).
     * This is a pure virtual function and hence MUST be
     * implemented by derived classes.
     */
    void init(amrex::AmrLevel &old) override;
    /**
     * Init data on this level after regridding if old AmrLevel
     * did not previously exist. This is a pure virtual function
     * and hence MUST be implemented by derived classes.
     */
    void init() override;
    /**
     * Do error estimation/tagging for regridding
     * Most examples should not need to override this and instead override
     * pre_tag_cells() and tag_cells()
     */
    void errorEst(amrex::TagBoxArray &a_tag_box_array, int a_clearval,
                  int a_tagval, amrex::Real a_time, int a_n_error_buf = 0,
                  int a_ngrow = 0) override;

    /**
     * Do any necessary work before tagging cells (e.g. calling FillPatch for
     * any variables for which derivatives are calculated).
     */
    virtual void pre_tag_cells() {}

    /**
     * Tag cells for regridding. This is a pure virtual function and hence MUST
     * be implemented by derived classes.
     */
    virtual void tag_cells(amrex::TagBoxArray &a_tag_box_array,
                           amrex::Real a_regrid_threshold) = 0;

    //! Do pre-plotfile work
    void writePlotFilePre(const std::string &dir,
                          std::ostream & /*os*/) override;

    //! Do post-plotfile work
    void writePlotFilePost(const std::string &dir,
                           std::ostream & /*os*/) override;

    //! Do pre-checkpoint work
    void checkPointPre(const std::string &a_dir, std::ostream &a_os) override;

    //! Do post-checkpoint work
    void checkPointPost(const std::string &a_dir, std::ostream &a_os) override;

    /// Virtual function for the problem specific parts of Advance
    virtual void specificAdvance() {}

    /// Virtual function for the problem specific parts of postTimeStep
    virtual void specificPostTimeStep() {}

    virtual void specificEvalRHS(amrex::MultiFab &a_soln,
                                 amrex::MultiFab &a_rhs,
                                 const amrex::Real a_time) = 0;

    virtual void specificUpdateODE(amrex::MultiFab & /*a_soln*/) {}

    //! Problem specific post restart
    virtual void specific_post_restart() {}

    //! Problem specific post init
    virtual void specific_post_init() {}

    //! Problem specific post-regrid
    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
    virtual void specific_post_regrid(int a_lbase, int a_new_finest) {}

    //! Problem specific pre plotfile
    virtual void specific_pre_plotfile(const std::string &a_dir,
                                       std::ostream &a_os)
    {
    }

    //! Problem specific post plotfile
    virtual void specific_post_plotfile(const std::string &a_dir,
                                        std::ostream &a_os)
    {
    }

    //! Problem specific pre checkpoint
    virtual void specific_pre_checkpoint(const std::string &a_dir,
                                         std::ostream &a_os)
    {
    }

    //! Problem specific post checkpoint
    virtual void specific_post_checkpoint(const std::string &a_dir,
                                          std::ostream &a_os)
    {
    }

    /// Returns true if m_time is the same as the time at the end of the current
    /// timestep on level a_level and false otherwise
    /// Useful to check whether to calculate something in postTimeStep (which
    /// might only be needed at the end of a_level's timestep)
    bool at_level_timestep_multiple(int a_level);

    BoundaryConditions m_boundaries; // the class for implementing BCs

    bool nan_check{};

  private:

    GRAMR *m_gramr_ptr = nullptr;
};

#endif /* GRAMRLEVEL_HPP_ */
