/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// Other includes
#include "BoundaryConditions.hpp"
#include "StateVariablesParmParse.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <numeric>
#include <string>

BoundaryConditions::params_t::params_t()
{
    // set defaults
    hi_boundary.fill(STATIC_BC);
    lo_boundary.fill(STATIC_BC);
    is_periodic.fill(true);

    vars_asymptotic_values.fill(0.0);
}

void BoundaryConditions::params_t::set_is_periodic(
    const std::array<int, AMREX_SPACEDIM> &a_is_periodic_int)
{
    FOR (idir)
    {
        is_periodic[idir] = static_cast<bool>(a_is_periodic_int[idir]);
        if (!is_periodic[idir])
        {
            nonperiodic_boundaries_exist = true;
        }
    }
}
void BoundaryConditions::params_t::set_hi_boundary(
    const std::array<int, AMREX_SPACEDIM> &a_hi_boundary)
{
    FOR (idir)
    {
        if (!is_periodic[idir])
        {
            hi_boundary[idir] = a_hi_boundary[idir];
            if (hi_boundary[idir] == REFLECTIVE_BC)
            {
                boundary_rhs_enforced       = true;
                boundary_solution_enforced  = true;
                reflective_boundaries_exist = true;
            }
            else if (hi_boundary[idir] == SOMMERFELD_BC)
            {
                boundary_rhs_enforced       = true;
                sommerfeld_boundaries_exist = true;
            }
            else if (hi_boundary[idir] == EXTRAPOLATING_BC)
            {
                boundary_rhs_enforced          = true;
                boundary_solution_enforced     = true;
                extrapolating_boundaries_exist = true;
            }
            else if (hi_boundary[idir] == MIXED_BC)
            {
                boundary_rhs_enforced          = true;
                boundary_solution_enforced     = true;
                sommerfeld_boundaries_exist    = true;
                extrapolating_boundaries_exist = true;
                mixed_boundaries_exist         = true;
            }
        }
    }
}
void BoundaryConditions::params_t::set_lo_boundary(
    const std::array<int, AMREX_SPACEDIM> &a_lo_boundary)
{
    FOR (idir)
    {
        if (!is_periodic[idir])
        {
            lo_boundary[idir] = a_lo_boundary[idir];
            if (lo_boundary[idir] == REFLECTIVE_BC)
            {
                boundary_solution_enforced  = true;
                reflective_boundaries_exist = true;
            }
            else if (lo_boundary[idir] == SOMMERFELD_BC)
            {
                boundary_rhs_enforced       = true;
                sommerfeld_boundaries_exist = true;
            }
            else if (lo_boundary[idir] == EXTRAPOLATING_BC)
            {
                boundary_rhs_enforced          = true;
                boundary_solution_enforced     = true;
                extrapolating_boundaries_exist = true;
            }
            else if (lo_boundary[idir] == MIXED_BC)
            {
                boundary_rhs_enforced          = true;
                boundary_solution_enforced     = true;
                sommerfeld_boundaries_exist    = true;
                extrapolating_boundaries_exist = true;
                mixed_boundaries_exist         = true;
            }
        }
    }
}

void BoundaryConditions::params_t::read_params(GRParmParse &pp)
{
    // still load even if not contained, to ensure printout saying parameters
    // were set to their default values
    std::array<int, AMREX_SPACEDIM> is_periodic_int{AMREX_D_DECL(1, 1, 1)};
    pp.load("isPeriodic", is_periodic_int, is_periodic_int);
    if (pp.contains("isPeriodic"))
    {
        set_is_periodic(is_periodic_int);
    }

    std::array<int, AMREX_SPACEDIM> hiBoundary{};
    pp.load("hi_boundary", hiBoundary, hi_boundary);
    if (pp.contains("hi_boundary"))
    {
        set_hi_boundary(hiBoundary);
    }

    std::array<int, AMREX_SPACEDIM> loBoundary{};
    pp.load("lo_boundary", loBoundary, lo_boundary);
    if (pp.contains("lo_boundary"))
    {
        set_lo_boundary(loBoundary);
    }

    if (sommerfeld_boundaries_exist)
    {
        size_t num_values = 0;
        std::vector<int> nonzero_asymptotic_vars;
        StateVariablesParmParse::load_vars_to_vector(
            pp, "nonzero_asymptotic_vars", nonzero_asymptotic_vars);
        const double default_value = 0.0;
        StateVariablesParmParse::load_values_to_array(
            pp, "nonzero_asymptotic_values", nonzero_asymptotic_vars,
            vars_asymptotic_values, default_value);
    }
    if (extrapolating_boundaries_exist)
    {
        pp.load("extrapolation_order", extrapolation_order, 1);
    }
    if (mixed_boundaries_exist)
    {
        std::vector<int> extrapolating_vars;
        StateVariablesParmParse::load_vars_to_vector(pp, "extrapolating_vars",
                                                     extrapolating_vars);
        for (int icomp = 0; icomp < NUM_VARS; icomp++)
        {
            bool is_extrapolating = false;
            // if the variable is not in extrapolating vars, it
            // is assumed to be sommerfeld by default
            for (auto &extrapolating_var : extrapolating_vars)
            {
                if (icomp == extrapolating_var)
                {
                    mixed_bc_vars_map.insert(
                        std::make_pair(icomp, EXTRAPOLATING_BC));
                    is_extrapolating = true;
                }
            }
            if (!is_extrapolating)
            {
                mixed_bc_vars_map.insert(std::make_pair(icomp, SOMMERFELD_BC));
            }
        }
    }
    if (nonperiodic_boundaries_exist)
    {
        // write out boundary conditions where non periodic - useful for
        // debug
        write_boundary_conditions(*this);
    }
}

/// define function sets members and is_defined set to true
void BoundaryConditions::define(std::array<double, AMREX_SPACEDIM> a_center,
                                const params_t &a_params,
                                const amrex::Geometry &a_geom, int a_num_ghosts)
{
    m_num_ghosts = a_num_ghosts;
    m_params     = a_params;
    FOR (i)
    {
        m_center[i] = a_center[i];
    }
    m_geom     = a_geom;
    is_defined = true;
}

/// change the asymptotic values of the variables for the Sommerfeld BCs
/// this will allow them to evolve during a simulation if necessary
void BoundaryConditions::set_vars_asymptotic_values(
    std::array<double, NUM_VARS> &vars_asymptotic_values)
{
    m_params.vars_asymptotic_values = vars_asymptotic_values;
}

void BoundaryConditions::write_reflective_conditions(int idir)
{
    amrex::Print()
        << "The variables that are parity odd in this direction are : " << '\n';
    for (int icomp = 0; icomp < NUM_VARS; icomp++)
    {
        int parity = get_state_var_parity(icomp, idir);
        if (parity == -1)
        {
            amrex::Print() << StateVariables::names[icomp] << "    ";
        }
    }
}

void BoundaryConditions::write_sommerfeld_conditions(int /*idir*/,
                                                     const params_t &a_params)
{
    amrex::Print() << "The non zero asymptotic values of the variables "
                      "in this direction are : "
                   << '\n';
    for (int icomp = 0; icomp < NUM_VARS; icomp++)
    {
        if (a_params.vars_asymptotic_values[icomp] != 0)
        {
            amrex::Print() << StateVariables::names[icomp] << " = "
                           << a_params.vars_asymptotic_values[icomp] << "    ";
        }
    }
    // not done for diagnostics
}

void BoundaryConditions::write_mixed_conditions(int idir,
                                                const params_t &a_params)
{
    // check all the vars have been assigned a BC - this should always be the
    // case because of how the params are assigned
    AMREX_ASSERT(a_params.mixed_bc_vars_map.size() == NUM_VARS);

    // now do the write out
    amrex::Print()
        << "The variables that use extrapolating bcs in this direction are : "
        << '\n';
    for (int icomp = 0; icomp < NUM_VARS; icomp++)
    {
        if (a_params.mixed_bc_vars_map.at(icomp) == EXTRAPOLATING_BC)
        {
            amrex::Print() << StateVariables::names[icomp] << "    ";
        }
    }
    amrex::Print() << '\n';
    amrex::Print()
        << "The other variables all use Sommerfeld boundary conditions."
        << '\n';
    write_sommerfeld_conditions(idir, a_params);
}

/// write out boundary params (used during setup for debugging)
void BoundaryConditions::write_boundary_conditions(const params_t &a_params)
{
    amrex::Print() << "You are using non periodic boundary conditions." << '\n';
    amrex::Print() << "The boundary params chosen are:  " << '\n';
    amrex::Print() << "---------------------------------" << '\n';

    std::map<int, std::string> bc_names = {
        {STATIC_BC,        "Static"       },
        {SOMMERFELD_BC,    "Sommerfeld"   },
        {REFLECTIVE_BC,    "Reflective"   },
        {EXTRAPOLATING_BC, "Extrapolating"},
        {MIXED_BC,         "Mixed"        }
    };
    FOR (idir)
    {
        if (!a_params.is_periodic[idir])
        {
            amrex::Print() << "- " << bc_names[a_params.hi_boundary[idir]]
                           << " boundaries in direction high " << idir << '\n';
            // high directions
            if (a_params.hi_boundary[idir] == REFLECTIVE_BC)
            {
                write_reflective_conditions(idir);
            }
            else if (a_params.hi_boundary[idir] == SOMMERFELD_BC)
            {
                write_sommerfeld_conditions(idir, a_params);
            }
            else if (a_params.hi_boundary[idir] == MIXED_BC)
            {
                write_mixed_conditions(idir, a_params);
            }
            amrex::Print() << "\n" << '\n';

            // low directions
            amrex::Print() << "- " << bc_names[a_params.lo_boundary[idir]]
                           << " boundaries in direction low " << idir << '\n';
            if (a_params.lo_boundary[idir] == REFLECTIVE_BC)
            {
                write_reflective_conditions(idir);
            }
            else if (a_params.lo_boundary[idir] == SOMMERFELD_BC)
            {
                write_sommerfeld_conditions(idir, a_params);
            }
            else if (a_params.lo_boundary[idir] == MIXED_BC)
            {
                write_mixed_conditions(idir, a_params);
            }
            amrex::Print() << "\n" << '\n';
        }
    }
    amrex::Print() << "---------------------------------" << '\n';
}

/// The function which returns the parity of each of the vars in
/// StateVariables.hpp (It is only required for reflective boundary conditions.)
// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
int BoundaryConditions::get_state_var_parity(int a_comp, int a_dir)
{
    BCParity comp_parity = StateVariables::parities[a_comp];

    if (((a_dir == 0) &&
         (comp_parity == BCParity::odd_x || comp_parity == BCParity::odd_xy ||
          comp_parity == BCParity::odd_xz ||
          comp_parity == BCParity::odd_xyz)) ||
        ((a_dir == 1) &&
         (comp_parity == BCParity::odd_y || comp_parity == BCParity::odd_xy ||
          comp_parity == BCParity::odd_yz ||
          comp_parity == BCParity::odd_xyz)) ||
        ((a_dir == 2) &&
         (comp_parity == BCParity::odd_z || comp_parity == BCParity::odd_xz ||
          comp_parity == BCParity::odd_yz || comp_parity == BCParity::odd_xyz)))
    {
        return -1;
    }
    return 1;
}

/// Get the boundary condition for given face
int BoundaryConditions::get_boundary_condition(amrex::Orientation face) const
{
    return face.isLow() ? m_params.lo_boundary[face.coordDir()]
                        : m_params.hi_boundary[face.coordDir()];
}

void BoundaryConditions::apply_sommerfeld_boundaries(
    amrex::MultiFab &a_rhs, const amrex::MultiFab &a_soln) const
{
    if (!m_params.sommerfeld_boundaries_exist)
    {
        return;
    }

    amrex::Vector<amrex::Box> sommboxes;
    {
        amrex::Box domain = m_geom.Domain();
        for (int idim = AMREX_SPACEDIM - 1; idim >= 0; --idim)
        {
            if (!m_params.is_periodic[idim])
            {
                int bclo = get_boundary_condition(
                    amrex::Orientation(idim, amrex::Orientation::low));
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(bclo != MIXED_BC,
                                                 "xxxx mixed bc todo");
                if (bclo == SOMMERFELD_BC)
                {
                    const int len  = domain.length(idim);
                    amrex::Box box = domain;
                    box.growHi(idim, -(len - m_num_ghosts));
                    domain.growLo(idim, -m_num_ghosts);
                    sommboxes.push_back(box);
                }
                int bchi = get_boundary_condition(
                    amrex::Orientation(idim, amrex::Orientation::high));
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(bchi != MIXED_BC,
                                                 "xxxx mixed bc todo");
                if (bchi == SOMMERFELD_BC)
                {
                    const int len  = domain.length(idim);
                    amrex::Box box = domain;
                    box.growLo(idim, -(len - m_num_ghosts));
                    domain.growHi(idim, -m_num_ghosts);
                    sommboxes.push_back(box);
                }
            }
        }
    }

    AMREX_ASSERT(amrex::almostEqual(m_geom.CellSize(0), m_geom.CellSize(1)) &&
                 amrex::almostEqual(m_geom.CellSize(0), m_geom.CellSize(2)));
    const auto dx     = m_geom.CellSize(0);
    amrex::Box domain = m_geom.Domain();
    for (amrex::OrientationIter orit; orit.isValid(); ++orit)
    {
        amrex::Orientation face = orit();
        int bc_on_face          = get_boundary_condition(face);
        if (m_geom.isPeriodic(face.coordDir()) || bc_on_face == REFLECTIVE_BC)
        {                      // xxxxx todo: what about other BCs?
            domain.grow(face); // to use the central derivative stencil
        }
    }
    const auto domlo  = domain.smallEnd();
    const auto domhi  = domain.bigEnd();
    const auto center = m_center;

    if (m_asymptotic_values.empty())
    {
        m_asymptotic_values.resize(NUM_VARS);
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, m_params.vars_asymptotic_values.begin(),
            m_params.vars_asymptotic_values.end(), m_asymptotic_values.begin());
    }
    auto *asymptotic_values = m_asymptotic_values.data();

#if defined(AMREX_USE_OMP) && !defined(AMREX_USE_GPU)
#pragma omp parallel
#endif
    for (amrex::MFIter mfi(a_rhs); mfi.isValid(); ++mfi)
    {
        const amrex::Box &valid_box                 = mfi.validbox();
        const amrex::Array4<amrex::Real const> &sol = a_soln.const_array(mfi);
        const amrex::Array4<amrex::Real> &rhs       = a_rhs.array(mfi);
        for (const auto &sommbox : sommboxes)
        {
            amrex::Box valid_sommbox = sommbox & valid_box;
            if (valid_sommbox.ok())
            {
                amrex::ParallelFor(
                    valid_sommbox, a_rhs.nComp(),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept
                    {
                        amrex::RealVect loc((i + 0.5) * dx - center[0],
                                            (j + 0.5) * dx - center[1],
                                            (k + 0.5) * dx - center[2]);
                        amrex::Real tmp = 0.;
                        amrex::IntVect iv(i, j, k);
                        for (int idir2 = 0; idir2 < AMREX_SPACEDIM; ++idir2)
                        {
                            amrex::IntVect iv_offset1 = iv;
                            amrex::IntVect iv_offset2 = iv;
                            amrex::Real d1            = NAN;
                            if (iv[idir2] == domlo[idir2])
                            {
                                iv_offset1[idir2] += +1;
                                iv_offset2[idir2] += +2;
                                d1 = (1.0 / dx) * (-1.5 * sol(iv, n) +
                                                   2.0 * sol(iv_offset1, n) -
                                                   0.5 * sol(iv_offset2, n));
                            }
                            else if (iv[idir2] == domhi[idir2])
                            {
                                iv_offset1[idir2] += -1;
                                iv_offset2[idir2] += -2;
                                d1 = (1.0 / dx) * (+1.5 * sol(iv, n) -
                                                   2.0 * sol(iv_offset1, n) +
                                                   0.5 * sol(iv_offset2, n));
                            }
                            else
                            {
                                iv_offset1[idir2] += +1;
                                iv_offset2[idir2] += -1;
                                d1                 = (0.5 / dx) *
                                     (sol(iv_offset1, n) - sol(iv_offset2, n));
                            }
                            // for each direction add dphidx * x^i
                            tmp += -d1 * loc[idir2];
                        }
                        // asymptotic values - these need to have been set in
                        // the params file
                        double radius =
                            std::sqrt(loc[0] * loc[0] + loc[1] * loc[1] +
                                      loc[2] * loc[2]);
                        rhs(i, j, k, n) =
                            (asymptotic_values[n] - sol(i, j, k, n) + tmp) *
                            (1. / radius);
                    });
            }
        }
    }
}

#if 0
//xxxxx
/// Fill the rhs boundary values appropriately based on the params set
void BoundaryConditions::fill_rhs_boundaries(const Side::LoHiSide a_side,
                                             const GRLevelData &a_soln,
                                             GRLevelData &a_rhs)
{
    AMREX_ASSERT(is_defined);
    BL_PROFILE("BoundaryConditions::fill_rhs_boundaries");

    // cycle through the directions, filling the rhs
    FOR(idir)
    {
        // only do something if this direction is not periodic
        if (!m_params.is_periodic[idir])
        {
            int boundary_condition = get_boundary_condition(a_side, idir);
            constexpr bool filling_rhs = true;
            fill_boundary_cells_dir(a_side, a_soln, a_rhs, idir,
                                    boundary_condition,
                                    Interval(0, NUM_VARS - 1),
                                    VariableType::state, filling_rhs);
        }
    }
}
#endif

#if 0
// xxxxx
/// fill solution boundary conditions, e.g. after interpolation
void BoundaryConditions::fill_solution_boundaries(const Side::LoHiSide a_side,
                                                  GRLevelData &a_state,
                                                  const Interval &a_comps)
{
    AMREX_ASSERT(is_defined);
    BL_PROFILE("BoundaryConditions::fill_solution_boundaries");

    // cycle through the directions
    FOR(idir)
    {
        // only do something if this direction is not periodic and solution
        // boundary enforced in this direction
        if (!m_params.is_periodic[idir])
        {
            int boundary_condition = get_boundary_condition(a_side, idir);

            // same copying of cells which we require for the rhs solution
            // but tell it we are not filling the rhs for mixed condition
            if ((boundary_condition == REFLECTIVE_BC) ||
                (boundary_condition == EXTRAPOLATING_BC) ||
                (boundary_condition == MIXED_BC))
            {
                const bool filling_rhs = false;
                fill_boundary_cells_dir(a_side, a_state, a_state, idir,
                                        boundary_condition, a_comps,
                                        VariableType::state, filling_rhs);
            }
        }
    }
}
#endif

#if 0
//xxxxx
/// fill diagnostic boundaries
void BoundaryConditions::fill_diagnostic_boundaries(const Side::LoHiSide a_side,
                                                    GRLevelData &a_state,
                                                    const Interval &a_comps)
{
    AMREX_ASSERT(is_defined);
    BL_PROFILE("BoundaryConditions::fill_diagnostic_boundaries");

    // cycle through the directions
    FOR(idir)
    {
        // only do something if this direction is not periodic
        if (!m_params.is_periodic[idir])
        {
            int boundary_condition = get_boundary_condition(a_side, idir);
            // for any non reflective BC, we just want to fill the ghosts with
            // something non nan so set the boundary condition to be
            // EXTRAPOLATING
            if (boundary_condition != REFLECTIVE_BC)
            {
                boundary_condition = EXTRAPOLATING_BC;
            }
            const bool filling_rhs = false;
            fill_boundary_cells_dir(a_side, a_state, a_state, idir,
                                    boundary_condition, a_comps,
                                    VariableType::derived, filling_rhs);
        }
    }
}
#endif

#if 0
//xxxxxx
/// Fill the boundary values appropriately based on the params set
/// in the direction dir
void BoundaryConditions::fill_boundary_cells_dir(
    const Side::LoHiSide a_side, const GRLevelData &a_soln, GRLevelData &a_out,
    const int dir, const int boundary_condition, const Interval &a_comps,
    const VariableType var_type, const bool filling_rhs)
{
    std::vector<int> comps_vector, sommerfeld_comps_vector,
        extrapolating_comps_vector;
    if (boundary_condition != MIXED_BC)
    {
        comps_vector.resize(a_comps.size());
        std::iota(comps_vector.begin(), comps_vector.end(), a_comps.begin());
    }
    else
    {
        for (int icomp = a_comps.begin(); icomp <= a_comps.end(); ++icomp)
        {
            if (m_params.mixed_bc_vars_map[icomp] == SOMMERFELD_BC)
                sommerfeld_comps_vector.push_back(icomp);
            else if (m_params.mixed_bc_vars_map[icomp] == EXTRAPOLATING_BC)
                extrapolating_comps_vector.push_back(icomp);
        }
    }

    // iterate through the boxes, shared amongst threads
    DataIterator dit = a_out.dataIterator();
    int nbox = dit.size();
#pragma omp parallel for default(shared)
    for (int ibox = 0; ibox < nbox; ++ibox)
    {
        DataIndex dind = dit[ibox];
        amrex::FArrayBox &out_box = a_out[dind];
        const amrex::FArrayBox &soln_box = a_soln[dind];
        amrex::Box this_box = out_box.box();
        amrex::IntVect offset_lo = -this_box.smallEnd() + m_domain_box.smallEnd();
        amrex::IntVect offset_hi = +this_box.bigEnd() - m_domain_box.bigEnd();

        // reduce box to the intersection of the box and the
        // problem domain ie remove all outer ghost cells
        this_box &= m_domain_box;
        // get the boundary box (may be Empty)
        amrex::Box boundary_box =
            get_boundary_box(a_side, dir, offset_lo, offset_hi, this_box);

        // now we have the appropriate box, fill it!
        BoxIterator bit(boundary_box);
        for (bit.begin(); bit.ok(); ++bit)
        {
            amrex::IntVect iv = bit();
            switch (boundary_condition)
            {
            // simplest case - boundary values are set to zero
            case STATIC_BC:
            {
                for (int icomp = a_comps.begin(); icomp <= a_comps.end();
                     ++icomp)
                {
                    out_box(iv, icomp) = 0.0;
                }
                break;
            }
            // Sommerfeld is outgoing radiation - only applies to rhs
            case SOMMERFELD_BC:
            {
                fill_sommerfeld_cell(out_box, soln_box, iv, comps_vector);
                break;
            }
            // Enforce a reflective symmetry in some direction
            case REFLECTIVE_BC:
            {
                fill_reflective_cell(out_box, iv, a_side, dir, comps_vector,
                                     var_type);
                break;
            }
            case EXTRAPOLATING_BC:
            {
                fill_extrapolating_cell(out_box, iv, a_side, dir, comps_vector,
                                        m_params.extrapolation_order);
                break;
            }
            case MIXED_BC:
            {
                fill_extrapolating_cell(out_box, iv, a_side, dir,
                                        extrapolating_comps_vector,
                                        m_params.extrapolation_order);
                if (filling_rhs)
                {
                    fill_sommerfeld_cell(out_box, soln_box, iv,
                                         sommerfeld_comps_vector);
                }
                break;
            }
            default:
                amrex::Abort(
                    "BoundaryCondition::Supplied boundary not supported.");
            } // end switch
        }     // end iterate over box
    }         // end iterate over boxes
}
#endif

void BoundaryConditions::fill_sommerfeld_cell(
    amrex::FArrayBox &rhs_box, const amrex::FArrayBox &soln_box,
    const amrex::IntVect a_iv, const std::vector<int> &sommerfeld_comps)
{
    amrex::Abort("xxxxx todo BoundaryConditions::fill_sommerfeld_cell");
    amrex::ignore_unused(rhs_box, soln_box, a_iv, sommerfeld_comps);
#if 0
//xxxxx
    // assumes an asymptotic value + radial waves and permits them
    // to exit grid with minimal reflections
    // get real position on the grid
    amrex::RealVect loc(a_iv + 0.5 * RealVect::Unit);
    loc *= m_dx;
    loc -= m_center;
    double radius_squared = 0.0;
    FOR(i) { radius_squared += loc[i] * loc[i]; }
    double radius = sqrt(radius_squared);
    amrex::IntVect lo_local_offset = a_iv - soln_box.smallEnd();
    amrex::IntVect hi_local_offset = soln_box.bigEnd() - a_iv;

    // Apply Sommerfeld BCs to each variable in sommerfeld_comps
    for (int icomp : sommerfeld_comps)
    {
        rhs_box(a_iv, icomp) = 0.0;
        FOR(idir2)
        {
            amrex::IntVect iv_offset1 = a_iv;
            amrex::IntVect iv_offset2 = a_iv;
            double d1;
            // bit of work to get the right stencils for near
            // the edges of the domain, only using second order
            // stencils for now
            if (lo_local_offset[idir2] < 1)
            {
                // near lo end
                iv_offset1[idir2] += +1;
                iv_offset2[idir2] += +2;
                d1 = 1.0 / m_dx *
                     (-1.5 * soln_box(a_iv, icomp) +
                      2.0 * soln_box(iv_offset1, icomp) -
                      0.5 * soln_box(iv_offset2, icomp));
            }
            else if (hi_local_offset[idir2] < 1)
            {
                // near hi end
                iv_offset1[idir2] += -1;
                iv_offset2[idir2] += -2;
                d1 = 1.0 / m_dx *
                     (+1.5 * soln_box(a_iv, icomp) -
                      2.0 * soln_box(iv_offset1, icomp) +
                      0.5 * soln_box(iv_offset2, icomp));
            }
            else
            {
                // normal case
                iv_offset1[idir2] += +1;
                iv_offset2[idir2] += -1;
                d1 =
                    0.5 / m_dx *
                    (soln_box(iv_offset1, icomp) - soln_box(iv_offset2, icomp));
            }

            // for each direction add dphidx * x^i / r
            rhs_box(a_iv, icomp) += -d1 * loc[idir2] / radius;
        }

        // asymptotic values - these need to have been set in
        // the params file
        rhs_box(a_iv, icomp) +=
            (m_params.vars_asymptotic_values[icomp] - soln_box(a_iv, icomp)) /
            radius;
    }
#endif
}

#if 0
//xxxxx
void BoundaryConditions::fill_extrapolating_cell(
    amrex::FArrayBox &out_box, const amrex::IntVect iv, const Side::LoHiSide a_side,
    const int dir, const std::vector<int> &extrapolating_comps,
    const int order) const
{
    for (int icomp : extrapolating_comps)
    {
        // current radius
        double radius = Coordinates<double>::get_radius(
            iv, m_dx, {m_center[0], m_center[1], m_center[2]});

        // vector of 2 nearest values and radii within the grid
        std::array<double, 2> value_at_point;
        std::array<double, 2> r_at_point;
        // how many units are we from domain boundary?
        int units_from_edge = 0;
        if (a_side == Side::Hi)
        {
            // how many units are we from domain boundary?
            units_from_edge = iv[dir] - m_domain_box.bigEnd(dir);
            // vector of 2 nearest values and radii within the grid
            for (int i = 0; i < 2; i++)
            {
                amrex::IntVect iv_tmp = iv;
                iv_tmp[dir] += -units_from_edge - i;
                FOR(idir)
                {
                    if (iv_tmp[idir] > m_domain_box.bigEnd(idir))
                    {
                        iv_tmp[idir] = m_domain_box.bigEnd(idir);
                    }
                    else if (iv_tmp[idir] < m_domain_box.smallEnd(idir))
                    {
                        iv_tmp[idir] = m_domain_box.smallEnd(idir);
                    }
                }
                value_at_point[i] = out_box(iv_tmp, icomp);
                r_at_point[i] = Coordinates<double>::get_radius(
                    iv_tmp, m_dx, {m_center[0], m_center[1], m_center[2]});
            }
        }
        else // Lo side
        {
            // how many units are we from domain boundary?
            units_from_edge = -iv[dir] + m_domain_box.smallEnd(dir);
            // vector of 2 nearest values within the grid
            for (int i = 0; i < 2; i++)
            {
                amrex::IntVect iv_tmp = iv;
                iv_tmp[dir] += units_from_edge + i;
                FOR(idir)
                {
                    if (iv_tmp[idir] > m_domain_box.bigEnd(idir))
                    {
                        iv_tmp[idir] = m_domain_box.bigEnd(idir);
                    }
                    else if (iv_tmp[idir] < m_domain_box.smallEnd(idir))
                    {
                        iv_tmp[idir] = m_domain_box.smallEnd(idir);
                    }
                }
                value_at_point[i] = out_box(iv_tmp, icomp);
                r_at_point[i] = Coordinates<double>::get_radius(
                    iv_tmp, m_dx, {m_center[0], m_center[1], m_center[2]});
            }
        }

        // assume some radial dependence and fit it
        double analytic_change = 0.0;
        // comp = const
        if (order == 0)
        {
            analytic_change = 0.0;
        }
        // comp = B + A*r
        else if (order == 1)
        {
            double delta_r_in_domain = r_at_point[1] - r_at_point[0];
            double A =
                (value_at_point[1] - value_at_point[0]) / delta_r_in_domain;
            double delta_r_here = radius - r_at_point[0];
            analytic_change = A * delta_r_here;
        }
        // other orders not supported yet
        else
        {
            amrex::Abort("Order not supported for boundary extrapolation.");
        }

        // set the value here to the extrapolated value
        out_box(iv, icomp) = value_at_point[0] + analytic_change;
    }
}
#endif

#if 0
//xxxxx
/// Copy the boundary values from src to dest
/// NB only acts if same box layout of input and output data
void BoundaryConditions::copy_boundary_cells(const Side::LoHiSide a_side,
                                             const GRLevelData &a_src,
                                             GRLevelData &a_dest)
{
    BL_PROFILE("BoundaryConditions::copy_boundary_cells");

    AMREX_ASSERT(is_defined);
    AMREX_ASSERT(a_src.nComp() == NUM_VARS);
    if (a_src.boxLayout() == a_dest.boxLayout())
    {
        // cycle through the directions
        FOR(idir)
        {
            // only do something if this direction is not periodic
            if (!m_params.is_periodic[idir])
            {
                // iterate through the boxes, shared amongst threads
                DataIterator dit = a_dest.dataIterator();
                int nbox = dit.size();
#pragma omp parallel for default(shared)
                for (int ibox = 0; ibox < nbox; ++ibox)
                {
                    DataIndex dind = dit[ibox];
                    amrex::FArrayBox &m_dest_box = a_dest[dind];
                    amrex::Box this_box = m_dest_box.box();
                    amrex::IntVect offset_lo =
                        -this_box.smallEnd() + m_domain_box.smallEnd();
                    amrex::IntVect offset_hi =
                        +this_box.bigEnd() - m_domain_box.bigEnd();

                    // reduce box to the intersection of the box and the
                    // problem domain ie remove all outer ghost cells
                    this_box &= m_domain_box;

                    // get the boundary box (may be Empty)
                    amrex::Box boundary_box = get_boundary_box(a_side, idir, offset_lo,
                                                        offset_hi, this_box);

                    BoxIterator bit(boundary_box);
                    for (bit.begin(); bit.ok(); ++bit)
                    {
                        amrex::IntVect iv = bit();
                        for (int icomp = 0; icomp < NUM_VARS; ++icomp)
                        {
                            m_dest_box(iv, icomp) = a_src[dind](iv, icomp);
                        }
                    } // end iterate over box
                }     // end iterate over boxes
            }         // end if(not periodic)
        }             // end iterate over spacedims
    }                 // end test for same box layout
}
#endif

#if 0
//xxxxx
/// Fill the fine boundary values in a_state
/// Required for interpolating onto finer levels at boundaries
void BoundaryConditions::interp_boundaries(GRLevelData &a_fine_state,
                                           GRLevelData &a_coarse_state,
                                           const Side::LoHiSide a_side)
{
    AMREX_ASSERT(is_defined);
    AMREX_ASSERT(a_fine_state.nComp() == NUM_VARS);
    AMREX_ASSERT(a_coarse_state.nComp() == NUM_VARS);
    BL_PROFILE("BoundaryConditions::interp_boundaries");

    // cycle through the directions
    FOR(idir)
    {
        // only do something if this direction is not periodic
        if (!m_params.is_periodic[idir])
        {
            // Ref ratio is always two
            int ref_ratio = 2;

            // create a coarsened fine layout and copy the coarse data onto
            // it
            DisjointBoxLayout coarsened_layout;
            coarsen(coarsened_layout, a_fine_state.disjointBoxLayout(),
                    ref_ratio * IntVect::Unit);
            GRLevelData coarsened_fine;
            coarsened_fine.define(coarsened_layout, NUM_VARS,
                                  m_num_ghosts * IntVect::Unit);
            amrex::Box coarse_domain_box = coarsen(m_domain_box, ref_ratio);

            // trick the copyTo into thinking the boundary cells are within
            // the domain by growing the domain
            amrex::Box grown_domain_box = coarse_domain_box;
            grown_domain_box.grow(m_num_ghosts * IntVect::Unit);
            Copier boundary_copier;
            boundary_copier.ghostDefine(
                a_coarse_state.disjointBoxLayout(),
                coarsened_fine.disjointBoxLayout(), grown_domain_box,
                m_num_ghosts * IntVect::Unit, m_num_ghosts * IntVect::Unit);
            a_coarse_state.copyTo(a_coarse_state.interval(), coarsened_fine,
                                  coarsened_fine.interval(), boundary_copier);

            // iterate through the coarse boxes, shared amongst threads
            DataIterator dit = coarsened_layout.dataIterator();
            int nbox = dit.size();
#pragma omp parallel for default(shared)
            for (int ibox = 0; ibox < nbox; ++ibox)
            {
                DataIndex dind = dit[ibox];
                amrex::FArrayBox &m_fine_box = a_fine_state[dind];
                amrex::FArrayBox &m_coarse_box = coarsened_fine[dind];
                amrex::Box this_box = m_coarse_box.box();
                amrex::Box fine_box = m_fine_box.box();
                amrex::IntVect offset_lo =
                    -this_box.smallEnd() + coarse_domain_box.smallEnd();
                amrex::IntVect offset_hi =
                    +this_box.bigEnd() - coarse_domain_box.bigEnd();

                // reduce box to the intersection of the box and the
                // problem domain ie remove all outer ghost cells
                this_box &= coarse_domain_box;

                // get the boundary box - remove one cell as we only want 2
                // coarse cells filled in each direction, to fill the 3 fine
                // cells on the level above
                amrex::Box boundary_box = get_boundary_box(a_side, idir, offset_lo,
                                                    offset_hi, this_box, 1);

                // define standard stencil for interp where not near
                // boundaries in other dirs
                amrex::IntVect default_offset =
                    IntVect::Zero + sign(a_side) * 2 * BASISV(idir);
                FourthOrderInterpStencil default_stencil(default_offset,
                                                         ref_ratio);

                // now interp the box from coarse to fine
                BoxIterator bit(boundary_box);
                for (bit.begin(); bit.ok(); ++bit)
                {
                    amrex::IntVect iv = bit();
                    amrex::IntVect lo_local_offset = iv - m_coarse_box.smallEnd();
                    amrex::IntVect hi_local_offset = m_coarse_box.bigEnd() - iv;

                    // bit of work to get the right stencils for near the
                    // edges of the box
                    bool near_boundary = false;
                    amrex::IntVect local_boundary_offset = IntVect::Zero;
                    FOR(idir2)
                    {
                        if (idir2 == idir)
                        {
                            local_boundary_offset[idir2] =
                                default_offset[idir2];
                        }
                        else if ((idir2 != idir) &&
                                 (lo_local_offset[idir2] > 1) &&
                                 (hi_local_offset[idir2] > 1))
                        {
                            local_boundary_offset[idir2] = 0;
                        }
                        else if ((idir2 != idir) &&
                                 (lo_local_offset[idir2] == 1))
                        {
                            local_boundary_offset[idir2] = -2;
                            near_boundary = true;
                        }
                        else if ((idir2 != idir) &&
                                 (hi_local_offset[idir2] == 1))
                        {
                            local_boundary_offset[idir2] = +2;
                            near_boundary = true;
                        }
                        else
                        {
                            amrex::Abort(
                                "BoundaryConditions::define bad boxes");
                        }
                    }

                    // if not near the boundary use the default stencil,
                    // otherwise use the one calculated locally
                    if (!near_boundary)
                    {
                        default_stencil.fillFine(m_fine_box, m_coarse_box, iv);
                    }
                    else
                    {
                        FourthOrderInterpStencil local_stencil(
                            local_boundary_offset, ref_ratio);
                        local_stencil.fillFine(m_fine_box, m_coarse_box, iv);
                    }
                } // end loop box
            }     // end loop boxes
        }         // end if is_periodic
    }             // end loop idir
}
#endif

#if 0
//xxxxx
/// get the boundary box to fill if we are at a boundary
Box BoundaryConditions::get_boundary_box(
    const Side::LoHiSide a_side, const int a_dir, const amrex::IntVect &offset_lo,
    const amrex::IntVect &offset_hi, amrex::Box &this_ghostless_box, int shrink_for_coarse)
{
    // default constructor gives empty box
    amrex::Box boundary_box;

    // check if we are over the edges of the domain - are we a boundary box?
    // if so create the box of the cells we want to fill
    if (((a_side == Side::Hi) && (offset_hi[a_dir] > 0)) ||
        ((a_side == Side::Lo) && (offset_lo[a_dir] > 0)))
    {
        // Get just the boundary box to iterate over, m_num_ghosts ghost
        // cells unless we are filling the coarse cells in the interp case
        // where we want to fill only two coarse ghost cells (to cover 3
        // fine ones)
        if (a_side == Side::Lo)
        {
            boundary_box = adjCellLo(this_ghostless_box, a_dir,
                                     m_num_ghosts - shrink_for_coarse);
        }
        else
        {
            boundary_box = adjCellHi(this_ghostless_box, a_dir,
                                     m_num_ghosts - shrink_for_coarse);
        }

        // adjust for any offsets - catches the corners etc
        // but only want to fill them once, so y fills x, z fills y and x
        // etc. Required in periodic direction corners in cases where there
        // are mixed boundaries, (otherwise these corners are full of nans)
        FOR(idir)
        {
            if (offset_lo[idir] > 0) // this direction is a low end boundary
            {
                if ((idir < a_dir) || (m_params.is_periodic[idir]))
                {
                    // grow it to fill the corners
                    boundary_box.growLo(idir, m_num_ghosts - shrink_for_coarse);
                }
            }
            else // cut off end ghost cell
            {
                if (idir != a_dir)
                {
                    boundary_box.growLo(idir, -shrink_for_coarse);
                }
            }

            if (offset_hi[idir] > 0) // this direction is a high end
                                     // boundary
            {
                if ((idir < a_dir) || (m_params.is_periodic[idir]))
                {
                    // grow it to fill the corners
                    boundary_box.growHi(idir, m_num_ghosts - shrink_for_coarse);
                }
            }
            else // cut off end ghost cell
            {
                if (idir != a_dir)
                {
                    boundary_box.growHi(idir, -shrink_for_coarse);
                }
            }
        }
    }
    return boundary_box;
}
#endif
