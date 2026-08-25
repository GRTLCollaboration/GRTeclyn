/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef BOUNDARYCONDITIONS_HPP_
#define BOUNDARYCONDITIONS_HPP_

// Our includes
#include "BCParity.hpp"
#include "DimensionDefinitions.hpp"
#include "GRParmParse.hpp"
#include "StateVariables.hpp"

#include <AMReX_MultiFab.H>

/// Class which deals with the boundaries at the edge of the physical domain in
/// cases where they are not periodic. Currently the options are first-order
/// extrapolation, Sommerfeld (outgoing radiation) and reflective BCs. The
/// conditions can differ in the high and low directions.
/// In cases where different variables/boundaries are required, the user should
/// (usually) write their own conditions class which inherits from this one.
/// Boundary handling combines AMReX ghost-cell filling with explicit RHS
/// updates where required. FIRST_ORDER_EXTRAPOLATION_BC and REFLECTIVE_BC are
/// translated into AMReX BCType values in GRAMRLevel::stateVariableSetUp(),
/// then imposed by AMReX when FillPatch fills the physical ghost cells.
class BoundaryConditions
{
  public:
    /// Boundary condition identifiers
    enum
    {
        UNSET_BC, // sentinel used where physical boundaries do not apply
        FIRST_ORDER_EXTRAPOLATION_BC, // AMReX foextrap
        SOMMERFELD_BC,
        REFLECTIVE_BC
    };

    /// Structure containing the boundary condition params
    struct params_t
    {
        std::array<int, AMREX_SPACEDIM> hi_condition{};
        std::array<int, AMREX_SPACEDIM> lo_condition{};
        std::array<bool, AMREX_SPACEDIM> is_periodic{};

        std::array<double, NUM_VARS> vars_asymptotic_values =
            StateVariables::asymptotic_values;
        params_t() = default;
        void fill_params();
        static void check_params();
        static std::array<int, AMREX_SPACEDIM>
        read_conditions(GRParmParse &a_boundary_pp, const char *a_name);
        [[nodiscard]] bool
        boundary_exists(const int a_boundary_condition) const;
    };

  protected:
    // Member values
    int m_num_ghosts{};       // the number of ghosts (usually 3)
    params_t m_params;        // the boundary params
    amrex::RealVect m_center; // the position of the center of the grid
    amrex::Geometry m_geom;   // the problem domain (excludes boundary cells)
    mutable amrex::Gpu::DeviceVector<double> m_asymptotic_values{};

  public:
    /// Default constructor - need to call define afterwards
    BoundaryConditions()  = default;
    ~BoundaryConditions() = default;

    BoundaryConditions(const BoundaryConditions &)            = delete;
    BoundaryConditions(BoundaryConditions &&)                 = delete;
    BoundaryConditions &operator=(const BoundaryConditions &) = delete;
    BoundaryConditions &operator=(BoundaryConditions &&)      = delete;

    /// Set the geometry-dependent members and load the boundary parameters
    void define(const amrex::Geometry &a_geom);

    /// change the asymptotic values of the variables for the Sommerfeld BCs
    /// this will allow them to evolve during a simulation if necessary
    void set_vars_asymptotic_values(
        std::array<double, NUM_VARS> &vars_asymptotic_values);

    /// The function which returns the parity of each of the vars in
    /// StateVariables.hpp (It is only required for reflective boundary
    /// conditions.)
    static int get_state_var_parity(int a_comp, int a_dir);

    /// Get the boundary condition for given face
    [[nodiscard]] int get_boundary_condition(amrex::Orientation face) const;

    /// Return whether the direction is periodic
    [[nodiscard]] bool is_periodic(int a_dir) const;

    /// Apply the Sommerfeld condition explicitly to the RHS after AMReX has
    /// filled its physical ghost cells using first-order extrapolation
    void apply_sommerfeld_boundaries(amrex::MultiFab &a_rhs,
                                     const amrex::MultiFab &a_soln) const;
};

#endif /* BOUNDARYCONDITIONS_HPP_ */
