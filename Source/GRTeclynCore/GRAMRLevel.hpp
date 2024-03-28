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
#include "UserVariables.hpp" // need NUM_VARS

#include <AMReX_AmrLevel.H>

#include <fstream>
#include <limits>
#include <sys/time.h>

enum StateType
{
    State_Type = 0,
    NUM_STATE_TYPE
};

// NOLINTNEXTLINE(cppcoreguidelines-special-member-functions)
class GRAMRLevel : public amrex::AmrLevel
{
  public:
    static void variableSetUp();
    static void variableCleanUp();

    GRAMRLevel();

    GRAMRLevel(amrex::Amr &papa, int lev, const amrex::Geometry &geom,
               const amrex::BoxArray &box_array,
               const amrex::DistributionMapping &distribution_mapping,
               amrex::Real time);

    ~GRAMRLevel() override;

    static const SimulationParameters &simParams();

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
    virtual void post_restart() override;
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
     * \brief Error estimation for regridding. This is a pure virtual
     * function and hence MUST be implemented by derived classes.
     * virtual void errorEst (amrex::TagBoxArray& tb, int clearval, int tagval,
     *                        amrex::Real time, int n_error_buf = 0,
     *                        int ngrow = 0);
     */

    //! Do pre-plotfile work
    void writePlotFilePre(const std::string &dir,
                          std::ostream & /*os*/) override;

    //! Do post-plotfile work
    void writePlotFilePost(const std::string &dir,
                           std::ostream & /*os*/) override;

    //! Return a MultiFab containing the derived data for this level.
    std::unique_ptr<amrex::MultiFab>
    derive(const std::string &name, amrex::Real time, int ngrow) override;

    //! Fill mf starting with the dcomp'th component with the derived quantity.
    void derive(const std::string &name, amrex::Real time,
                amrex::MultiFab &multifab, int dcomp) override;

    /// Virtual function for the problem specific parts of Advance
    virtual void specificAdvance() {}

    /// Virtual function for the problem specific parts of postTimeStep
    virtual void specificPostTimeStep() {}

    virtual void specificEvalRHS(amrex::MultiFab &a_soln,
                                 amrex::MultiFab &a_rhs,
                                 const double a_time) = 0;

    virtual void specificUpdateODE(amrex::MultiFab & /*a_soln*/) {}

    BoundaryConditions m_boundaries; // the class for implementing BCs

    int m_verbosity = 0; //!< Level of verbosity of the output
    int m_num_ghosts{};  //!< Number of ghost cells
    bool m_is_writing_plotfile = false;

    // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
    static amrex::Vector<std::string> plot_constraints;

  private:

    GRAMR *m_gramr_ptr = nullptr;
};

#endif /* GRAMRLEVEL_HPP_ */
