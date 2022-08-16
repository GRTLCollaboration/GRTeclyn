/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef GRAMRLEVEL_HPP_
#define GRAMRLEVEL_HPP_

// Other includes
#include "BoundaryConditions.hpp"
#include "GRAMR.hpp"
#include "GRLevelData.hpp"
#include "InterpSource.hpp"
#include "SimulationParameters.hpp"
#include "UserVariables.hpp" // need NUM_VARS

// Chombo includes
// xxxxx#include "AMRLevel.H"
//#include "CoarseAverage.H"
//#include "FourthOrderFillPatch.H"
//#include "LevelFluxRegister.H" //We don't actually use flux conservation but
//Chombo assumes we do #include "LevelRK4.H" #include "LoadBalance.H"

#include <fstream>
#include <limits>
#include <sys/time.h>

#include <AMReX_AmrLevel.H>

enum StateType
{
    State_Type = 0,
    NUM_STATE_TYPE
};

class GRAMRLevel : public amrex::AmrLevel // xxxxx, public InterpSource
{
  public:
    static void variableSetUp();
    static void variableCleanUp();

    GRAMRLevel();

    GRAMRLevel(amrex::Amr &papa, int lev, const amrex::Geometry &geom,
               const amrex::BoxArray &ba, const amrex::DistributionMapping &dm,
               amrex::Real time);

    virtual ~GRAMRLevel();

    SimulationParameters const &simParams() const;

    /**
     * \brief Compute the initial time step.
     */
    virtual void computeInitialDt(
        int finest_level, int sub_cycle, amrex::Vector<int> &n_cycle,
        const amrex::Vector<amrex::IntVect> &ref_ratio,
        amrex::Vector<amrex::Real> &dt_level, amrex::Real stop_time) override;
    /**
     * \brief Compute the next time step.
     */
    virtual void computeNewDt(int finest_level, int sub_cycle,
                              amrex::Vector<int> &n_cycle,
                              const amrex::Vector<amrex::IntVect> &ref_ratio,
                              amrex::Vector<amrex::Real> &dt_min,
                              amrex::Vector<amrex::Real> &dt_level,
                              amrex::Real stop_time,
                              int post_regrid_flag) override;
    /**
     * \brief Do an integration step on this level.  Returns maximum safe
     * time step.  This is a pure virtual function and hence MUST
     * be implemented by derived classes.
     */
    virtual amrex::Real advance(amrex::Real time, amrex::Real dt, int iteration,
                                int ncycle) override;

    /**
     * \brief Contains operations to be done after a timestep.  This is a
     * pure virtual function and hence MUST be implemented by derived
     * classes.
     */
    virtual void post_timestep(int iteration) override;
    /**
     * \brief Operations to be done after regridding
     * This is a pure virtual function and hence MUST be
     * implemented by derived classes.
     */
    virtual void post_regrid(int lbase, int new_finest) override;
    /**
     * \brief Operations to be done after initialization.
     * This is a pure virtual function and hence MUST be
     * implemented by derived classes.
     */
    virtual void post_init(amrex::Real stop_time) override;
    /**
     * \brief Init data on this level from another AmrLevel (during regrid).
     * This is a pure virtual function and hence MUST be
     * implemented by derived classes.
     */
    virtual void init(amrex::AmrLevel &old) override;
    /**
     * Init data on this level after regridding if old AmrLevel
     * did not previously exist. This is a pure virtual function
     * and hence MUST be implemented by derived classes.
     */
    virtual void init() override;
    /**
     * \brief Error estimation for regridding. This is a pure virtual
     * function and hence MUST be implemented by derived classes.
     */
    virtual void errorEst(amrex::TagBoxArray &tb, int clearval, int tagval,
                          amrex::Real time, int n_error_buf = 0,
                          int ngrow = 0) override;

    // xxxxx  public:
    // xxxxx    GRAMRLevel(GRAMR &gr_amr, const SimulationParameters &a_p, int
    // a_verbosity); xxxxx xxxxx xxxxx  public: xxxxx    /// Do casting from
    // AMRLevel to GRAMRLevel and stop if this isn't possible xxxxx    static
    // const GRAMRLevel *gr_cast(const AMRLevel *const amr_level_ptr); xxxxx
    // static GRAMRLevel *gr_cast(AMRLevel *const amr_level_ptr); xxxxx xxxxx
    // const GRLevelData & xxxxx    getLevelData(const VariableType var_type =
    // VariableType::evolution) const; xxxxx xxxxx    bool contains(const
    // std::array<double, AMREX_SPACEDIM> &point) const; xxxxx xxxxx  private:
    // xxxxx    // define
    // xxxxx    virtual void define(AMRLevel *a_coarser_level_ptr,
    // xxxxx                        const Box &a_problem_domain, int a_level,
    // xxxxx                        int a_ref_ratio);
    // xxxxx
    // xxxxx    // define
    // xxxxx    virtual void define(AMRLevel *a_coarser_level_ptr,
    // xxxxx                        const ProblemDomain &a_problem_domain, int
    // a_level, xxxxx                        int a_ref_ratio); xxxxx xxxxx    ///
    // advance by one timestep xxxxx    virtual Real advance(); xxxxx xxxxx ///
    // things to do after a timestep xxxxx    virtual void postTimeStep(); xxxxx
    // xxxxx    /// things to do before tagging cells (e.g. filling ghosts)
    // xxxxx    virtual void preTagCells();
    // xxxxx
    // xxxxx    /// tag cells that need to be refined
    // xxxxx    virtual void tagCells(IntVectSet &a_tags);
    // xxxxx
    // xxxxx    /// create tags at initialization
    // xxxxx    virtual void tagCellsInit(IntVectSet &a_tags);
    // xxxxx
    // xxxxx    /// regrid
    // xxxxx    virtual void regrid(const Vector<Box> &a_new_grids);
    // xxxxx
    // xxxxx    /// things to do after regridding
    // xxxxx    virtual void postRegrid(int a_base_level);
    // xxxxx
    // xxxxx    /// initialize grids
    // xxxxx    virtual void initialGrid(const Vector<Box> &a_new_grids);
    // xxxxx
    // xxxxx    /// things to do after initialization
    // xxxxx    virtual void postInitialize();
    // xxxxx
    // xxxxx    /// compute the size of the timestep
    // xxxxx    virtual Real computeDt();
    // xxxxx
    // xxxxx    /// compute the size of the initial timestep
    // xxxxx    virtual Real computeInitialDt();
    // xxxxx
    // xxxxx    DisjointBoxLayout loadBalance(const Vector<Box> &a_grids);
    // xxxxx
    // xxxxx#ifdef AMREX_USE_HDF5
    // xxxxx    virtual void writeCheckpointHeader(HDF5Handle &a_handle) const;
    // xxxxx
    // xxxxx    virtual void writeCheckpointLevel(HDF5Handle &a_handle) const;
    // xxxxx
    // xxxxx    virtual void readCheckpointHeader(HDF5Handle &a_handle);
    // xxxxx
    // xxxxx    virtual void readCheckpointLevel(HDF5Handle &a_handle);
    // xxxxx
    // xxxxx    virtual void writePlotHeader(HDF5Handle &a_handle) const;
    // xxxxx
    // xxxxx    virtual void writePlotLevel(HDF5Handle &a_handle) const;
    // xxxxx#endif
    // xxxxx
    // xxxxx  public:
    // xxxxx    /// evaluate d(soln)/dt at current time based on soln
    // xxxxx    void evalRHS(GRLevelData &rhs,          //!< d(soln)/dt based on
    // soln xxxxx                 GRLevelData &soln,         //!< soln at current
    // time xxxxx                 LevelFluxRegister &fineFR, //!< flux register
    // w/ finer level xxxxx                 LevelFluxRegister &crseFR, //!< flux
    // register w/ crse level xxxxx                 const GRLevelData
    // &oldCrseSoln, //!< old-time crse solution xxxxx                 Real
    // oldCrseTime,               //!< old crse time xxxxx                 const
    // GRLevelData &newCrseSoln, //!< new-time crse solution xxxxx Real
    // newCrseTime,               //!< new crse time xxxxx                 Real
    // time,      //!< current time centering of soln xxxxx                 Real
    // fluxWeight //!< weight to apply to fluxRegister updates xxxxx    ); xxxxx
    // xxxxx    /// implements soln += dt*rhs
    // xxxxx    void updateODE(GRLevelData &soln, const GRLevelData &rhs, Real
    // dt); xxxxx xxxxx    /// define data holder newSoln based on existingSoln,
    // including ghost cell xxxxx    /// specification xxxxx    void
    // defineSolnData(GRLevelData &newSoln, const GRLevelData &existingSoln);
    // xxxxx
    // xxxxx    /// define data holder for RHS based on existingSoln including
    // ghost cell xxxxx    /// specification (which in most cases is no ghost
    // cells) xxxxx    void defineRHSData(GRLevelData &newRHS, const GRLevelData
    // &existingSoln); xxxxx xxxxx    /// copy data from src into dest xxxxx void
    // copySolnData(GRLevelData &dest, const GRLevelData &src); xxxxx xxxxx ///
    // Virtual function for the problem specific parts of Advance xxxxx virtual
    // void specificAdvance() {} xxxxx xxxxx    /// Virtual function for the
    // problem specific parts of postTimeStep xxxxx    virtual void
    // specificPostTimeStep() {} xxxxx xxxxx    /// (Pure) virtual function for
    // the initial data calculation xxxxx    virtual void initialData() = 0;
    // xxxxx
    // xxxxx    /// Computes which cells have insufficient resolution and should
    // be tagged xxxxx    virtual void computeTaggingCriterion(FArrayBox
    // &tagging_criterion, xxxxx                                         const
    // FArrayBox &current_state) = 0; xxxxx xxxxx#ifdef AMREX_USE_HDF5 xxxxx ///
    // Things to do immediately before checkpointing xxxxx    virtual void
    // preCheckpointLevel() {} xxxxx xxxxx    /// Things to do immediately before
    // writing plot files xxxxx    virtual void prePlotLevel() {} xxxxx xxxxx ///
    // Things to do immediately after restart from checkpoint xxxxx    virtual
    // void postRestart() {} xxxxx#endif xxxxx xxxxx    virtual void
    // specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs, xxxxx const
    // double a_time) = 0; xxxxx xxxxx    virtual void
    // specificUpdateODE(GRLevelData &a_soln, xxxxx const GRLevelData &a_rhs,
    // Real a_dt) xxxxx    { xxxxx    } xxxxx xxxxx    double get_dx() const;
    // xxxxx
    // xxxxx    /// Returns true if m_time is the same as the time at the end of
    // the current xxxxx    /// timestep on level a_level and false otherwise
    // xxxxx    /// Useful to check whether to calculate something in
    // postTimeStep (which xxxxx    /// might only be needed at the end of
    // a_level's timestep) xxxxx    bool at_level_timestep_multiple(int a_level)
    // const; xxxxx xxxxx    /// Fill all [either] evolution or diagnostic ghost
    // cells xxxxx    virtual void fillAllGhosts( xxxxx        const VariableType
    // var_type = VariableType::evolution, xxxxx        const Interval &a_comps =
    // Interval(0, std::numeric_limits<int>::max())); xxxxx xxxxx  protected:
    // xxxxx    /// Fill all evolution ghosts cells (i.e. those in m_state_new)
    // xxxxx    virtual void
    // xxxxx    fillAllEvolutionGhosts(const Interval &a_comps = Interval(0,
    // NUM_VARS - 1)); xxxxx xxxxx    /// Fill all diagnostics ghost cells (i.e.
    // those in m_state_diagnostics) xxxxx    virtual void
    // fillAllDiagnosticsGhosts( xxxxx        const Interval &a_comps =
    // Interval(0, NUM_DIAGNOSTIC_VARS - 1)); xxxxx xxxxx    /// Fill ghosts
    // cells from boxes on this level only. Do not interpolate xxxxx    ///
    // between levels. xxxxx    virtual void xxxxx    fillIntralevelGhosts(const
    // Interval &a_comps = Interval(0, NUM_VARS - 1)); xxxxx xxxxx    /// This
    // function is used to fill ghost cells outside the domain xxxxx    /// (for
    // non-periodic boundary conditions, where values depend on state) xxxxx
    // virtual void fillBdyGhosts(GRLevelData &a_state, xxxxx const Interval
    // &a_comps = Interval(0, NUM_VARS - xxxxx 1)); xxxxx xxxxx    /// This
    // function is used to copy ghost cells outside the domain xxxxx    /// (for
    // non-periodic boundary conditions, where boundaries evolve via rhs) xxxxx
    // virtual void copyBdyGhosts(const GRLevelData &a_src, GRLevelData &a_dest);
    // xxxxx
    // xxxxx    /// This function is used to define the exchange copiers
    // required for xxxxx    /// copying ghost cells between boxes xxxxx virtual
    // void defineExchangeCopier(const DisjointBoxLayout &a_level_domain); xxxxx
    // xxxxx    void printProgress(const std::string &from) const;
    // xxxxx
    // xxxxx    BoundaryConditions m_boundaries; // the class for implementing
    // BCs xxxxx xxxxx    GRLevelData m_state_old; //!< the solution at the old
    // time xxxxx    GRLevelData m_state_new; //!< the solution at the new time
    // xxxxx    GRLevelData m_state_diagnostics;
    // xxxxx    Real m_dx; //!< grid spacing
    // xxxxx    double m_restart_time;
    // xxxxx
    // xxxxx    GRAMR &m_gr_amr; //!< The GRAMR object containing this
    // GRAMRLevel xxxxx xxxxx    // params xxxxx    SimulationParameters m_p;
    // //!< Holds parameters necessary for the simulation
    int m_verbosity; //!< Level of verbosity of the output
    // xxxxx
    // xxxxx    Copier m_exchange_copier; //!< copier (for ghost cells on same
    // level) xxxxx xxxxx    CoarseAverage m_coarse_average; //!< Averages from
    // fine to coarse level xxxxx xxxxx    FourthOrderFillPatch m_patcher; //!<
    // Organises interpolation from coarse to xxxxx //!< fine levels of ghosts
    // xxxxx    FourthOrderFillPatch
    // xxxxx        m_patcher_diagnostics; //!< Organises interpolation from
    // coarse to xxxxx                               //!< fine levels of ghosts
    // for diagnostics xxxxx    FourthOrderFineInterp m_fine_interp; //!<
    // executes the interpolation from xxxxx //!< coarse to fine when regridding
    // xxxxx
    // xxxxx    DisjointBoxLayout m_grids;       //!< Holds grid setup (the
    // layout of boxes) xxxxx    DisjointBoxLayout m_grown_grids; //!< Holds
    // grown grid setup (for xxxxx                                     //!<
    // Sommerfeld BCs) xxxxx xxxxx  public: xxxxx    const int m_num_ghosts; //!<
    // Number of ghost cells
};

#endif /* GRAMRLEVEL_HPP_ */
