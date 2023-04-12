/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BOUNDARYCONDITIONS_HPP_
#define BOUNDARYCONDITIONS_HPP_

// Our includes
#include "DimensionDefinitions.hpp"
#include "GRParmParse.hpp"
#include "UserVariables.hpp"
#include "VariableType.hpp"

#include <AMReX_MultiFab.H>

/// Class which deals with the boundaries at the edge of the physical domain in
/// cases where they are not periodic. Currently only options are static BCs,
/// sommerfeld (outgoing radiation) and reflective. The conditions can differ in
/// the high and low directions.
/// In cases where different variables/boundaries are required, the user should
/// (usually) write their own conditions class which inherits from this one.
/// Note that these conditions enforce a certain rhs based on the current values
/// of the grid variables. (Another option would be to enforce grid values, e.g.
/// by extrapolating from within the grid.)
class BoundaryConditions
{
  public:
    /// enum for possible boundary states
    enum
    {
        STATIC_BC,
        SOMMERFELD_BC,
        REFLECTIVE_BC,
        EXTRAPOLATING_BC,
        MIXED_BC
    };

    /// enum for possible parity states
    enum
    {
        EVEN,
        ODD_X,
        ODD_Y,
        ODD_Z,
        ODD_XY,
        ODD_YZ,
        ODD_XZ,
        ODD_XYZ,
        UNDEFINED
    };

    /// Structure containing the boundary condition params
    struct params_t
    {
        std::array<int, AMREX_SPACEDIM> hi_boundary{};
        std::array<int, AMREX_SPACEDIM> lo_boundary{};
        std::array<bool, AMREX_SPACEDIM> is_periodic{};
        bool nonperiodic_boundaries_exist{false};
        bool boundary_solution_enforced{false};
        bool boundary_rhs_enforced{false};
        bool reflective_boundaries_exist{false};
        bool sommerfeld_boundaries_exist{false};
        bool extrapolating_boundaries_exist{};
        bool mixed_boundaries_exist{};

        std::array<int, NUM_VARS> vars_parity{};
        std::array<int, NUM_DIAGNOSTIC_VARS>
            vars_parity_diagnostic{}; /* needed only in AMRInterpolator */
        std::array<double, NUM_VARS> vars_asymptotic_values{};
        std::map<int, int> mixed_bc_vars_map;
        int extrapolation_order{1};
        params_t(); // sets the defaults
        void
        set_is_periodic(const std::array<bool, AMREX_SPACEDIM> &a_is_periodic);
        void
        set_hi_boundary(const std::array<int, AMREX_SPACEDIM> &a_hi_boundary);
        void
        set_lo_boundary(const std::array<int, AMREX_SPACEDIM> &a_lo_boundary);
        void read_params(GRParmParse &pp);
    };

  protected:
    // Member values
    int m_num_ghosts{};       // the number of ghosts (usually 3)
    params_t m_params;        // the boundary params
    amrex::RealVect m_center; // the position of the center of the grid
    amrex::Geometry m_geom;   // the problem domain (excludes boundary cells)
    bool is_defined{
        false}; // whether the BoundaryConditions class members are defined
    mutable amrex::Gpu::DeviceVector<double> m_asymptotic_values;

  public:
    /// Default constructor - need to call define afterwards
    BoundaryConditions()  = default;
    ~BoundaryConditions() = default;

    BoundaryConditions(const BoundaryConditions &)            = delete;
    BoundaryConditions(BoundaryConditions &&)                 = delete;
    BoundaryConditions &operator=(const BoundaryConditions &) = delete;
    BoundaryConditions &operator=(BoundaryConditions &&)      = delete;

    /// define function sets members and is_defined set to true
    void define(std::array<double, AMREX_SPACEDIM> a_center,
                const params_t &a_params, const amrex::Geometry &a_geom,
                int a_num_ghosts);

    /// change the asymptotic values of the variables for the Sommerfeld BCs
    /// this will allow them to evolve during a simulation if necessary
    void set_vars_asymptotic_values(
        std::array<double, NUM_VARS> &vars_asymptotic_values);

    /// write out boundary params (used during setup for debugging)
    static void write_boundary_conditions(const params_t &a_params);

    /// The function which returns the parity of each of the vars in
    /// UserVariables.hpp The parity should be defined in the params file, and
    /// will be output to the pout files for checking at start/restart of
    /// simulation (It is only required for reflective boundary conditions.)
    int
    get_var_parity(int a_comp, int a_dir,
                   const VariableType var_type = VariableType::evolution) const;

    /// static version used for initial output of boundary values
    static int
    get_var_parity(int a_comp, int a_dir, const params_t &a_params,
                   const VariableType var_type = VariableType::evolution);

    /// Get the boundary condition for given face
    int get_boundary_condition(amrex::Orientation face) const;

    /// Apply Sommerfeld BC to RHS
    void apply_sommerfeld_boundaries(amrex::MultiFab &a_rhs,
                                     const amrex::MultiFab &a_soln) const;

#if 0
//xxxxx
    /// Fill the rhs boundary values appropriately based on the params set
    void fill_rhs_boundaries(const Side::LoHiSide a_side,
                             const GRLevelData &a_soln, GRLevelData &a_rhs);

    /// enforce solution boundary conditions, e.g. after interpolation
    void fill_solution_boundaries(
        const Side::LoHiSide a_side, GRLevelData &a_state,
        const Interval &a_comps = Interval(0, NUM_VARS - 1));

    /// fill diagnostic boundaries - used in AMRInterpolator
    void fill_diagnostic_boundaries(
        const Side::LoHiSide a_side, GRLevelData &a_state,
        const Interval &a_comps = Interval(0, NUM_DIAGNOSTIC_VARS - 1));

    /// Fill the boundary values appropriately based on the params set
    /// in the direction dir
    void fill_boundary_cells_dir(const Side::LoHiSide a_side,
                                 const GRLevelData &a_soln, GRLevelData &a_out,
                                 const int dir, const int boundary_condition,
                                 const Interval &a_comps,
                                 const VariableType var_type,
                                 const bool filling_rhs);

    /// Copy the boundary values from src to dest
    /// NB assumes same box layout of input and output data
    void copy_boundary_cells(const Side::LoHiSide a_side,
                             const GRLevelData &a_src, GRLevelData &a_dest);

    /// Fill the fine boundary values in a_state
    /// Required for interpolating onto finer levels at boundaries
    void interp_boundaries(GRLevelData &a_fine_state,
                           GRLevelData &a_coarse_state,
                           const Side::LoHiSide a_side);

    /// Get the boundary condition for a_dir and a_side
    int get_boundary_condition(const Side::LoHiSide a_side, const int a_dir);

    /// get the boundary box to fill if we are at a boundary
    amrex::Box get_boundary_box(const Side::LoHiSide a_side, const int a_dir,
                         const amrex::IntVect &offset_lo, const amrex::IntVect &offset_hi,
                         amrex::Box &this_ghostless_box, int shrink_for_coarse = 0);

    /// This function takes a default constructed open DisjointBoxLayout and
    /// grows the boxes lying along the boundary to include the boundaries if
    /// necessary (i.e. in the Sommerfeld BC case). It is used to define the
    /// correct DisjointBoxLayout for the exchange copier so that shared
    /// boundary ghosts are exchanged correctly.
//xxxxx    void expand_grids_to_boundaries(DisjointBoxLayout &a_out_grids,
//                                    const DisjointBoxLayout &a_in_grids);
#endif

    friend class ExpandGridsToBoundaries;

  private:
    /// write out reflective conditions
    static void write_reflective_conditions(int idir, const params_t &a_params);

    /// write out sommerfeld conditions
    static void write_sommerfeld_conditions(int idir, const params_t &a_params);

    /// write out mixed conditions
    static void write_mixed_conditions(int idir, const params_t &a_params);

    static void fill_sommerfeld_cell(amrex::FArrayBox &rhs_box,
                                     const amrex::FArrayBox &soln_box,
                                     const amrex::IntVect a_iv,
                                     const std::vector<int> &sommerfeld_comps);

#if 0
//xxxxx    void fill_extrapolating_cell(amrex::FArrayBox &out_box, const amrex::IntVect iv,
                                 const Side::LoHiSide a_side, const int dir,
                                 const std::vector<int> &extrapolating_comps,
                                 const int order = 1) const;

    void fill_reflective_cell(
        amrex::FArrayBox &out_box, const amrex::IntVect iv, const Side::LoHiSide a_side,
        const int dir, const std::vector<int> &reflective_comps,
        const VariableType var_type = VariableType::evolution) const;
#endif
};

#endif /* BOUNDARYCONDITIONS_HPP_ */
