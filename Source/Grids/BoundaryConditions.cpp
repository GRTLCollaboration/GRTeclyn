/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// Other includes
#include "BoundaryConditions.hpp"

#include <array>
#include <cmath>
#include <map>
#include <string>

namespace
{
// NOLINTNEXTLINE(bugprone-throwing-static-initialization)
const std::map<std::string, int> boundary_conditions_by_name = {
    {"UNSET_BC",                     BoundaryConditions::UNSET_BC     },
    {"FIRST_ORDER_EXTRAPOLATION_BC",
     BoundaryConditions::FIRST_ORDER_EXTRAPOLATION_BC                 },
    {"SOMMERFELD_BC",                BoundaryConditions::SOMMERFELD_BC},
    {"REFLECTIVE_BC",                BoundaryConditions::REFLECTIVE_BC}
};

} // namespace

std::array<int, AMREX_SPACEDIM>
BoundaryConditions::params_t::read_conditions(GRParmParse &a_boundary_pp,
                                              const char *a_name)
{
    std::array<std::string, AMREX_SPACEDIM> condition_names{};
    a_boundary_pp.get(a_name, condition_names);

    std::array<int, AMREX_SPACEDIM> conditions{};
    FOR (idir)
    {
        const auto boundary_condition_entry =
            boundary_conditions_by_name.find(condition_names[idir]);
        if (boundary_condition_entry == boundary_conditions_by_name.end())
        {
            a_boundary_pp.error(
                a_name,
                "entries must be UNSET_BC, FIRST_ORDER_EXTRAPOLATION_BC, "
                "SOMMERFELD_BC, or REFLECTIVE_BC");
        }
        else
        {
            // The map value is the integer identifier for the named boundary
            // condition.
            conditions[idir] = boundary_condition_entry->second;
        }
    }
    return conditions;
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
bool BoundaryConditions::params_t::boundary_exists(
    const int a_boundary_condition) const
{
    FOR (idir)
    {
        if (!is_periodic[idir] && (hi_condition[idir] == a_boundary_condition ||
                                   lo_condition[idir] == a_boundary_condition))
        {
            return true;
        }
    }
    return false;
}

void BoundaryConditions::params_t::fill_params()
{
    // still load even if not contained, to ensure printout saying parameters
    // were set to their default values
    GRParmParse boundary_pp("boundary");
    GRParmParse geom_pp("geometry");
    GRParmParse pp;

    std::array<int, AMREX_SPACEDIM> is_periodic_int = {0, 0, 0};
    geom_pp.get("is_periodic", is_periodic_int);

    FOR (idir)
    {
        this->is_periodic[idir] = static_cast<bool>(is_periodic_int[idir]);
    }

    this->hi_condition = read_conditions(boundary_pp, "hi_condition");
    this->lo_condition = read_conditions(boundary_pp, "lo_condition");
}

void BoundaryConditions::params_t::check_params()
{
    // Register defaults so that the resolved values appear in the parameter
    // output even when they were not supplied by the user.
    GRParmParse boundary_pp("boundary");
    GRParmParse geom_pp("geometry");

    // set defaults
    std::array<std::string, AMREX_SPACEDIM> hi_condition_names{};
    hi_condition_names.fill("UNSET_BC");
    boundary_pp.queryAdd("hi_condition", hi_condition_names);
    const auto hi_conditions = read_conditions(boundary_pp, "hi_condition");

    std::array<std::string, AMREX_SPACEDIM> lo_condition_names{};
    lo_condition_names.fill("UNSET_BC");
    boundary_pp.queryAdd("lo_condition", lo_condition_names);
    const auto lo_conditions = read_conditions(boundary_pp, "lo_condition");

    std::array<int, AMREX_SPACEDIM> is_periodic_int{};
    geom_pp.get("is_periodic", is_periodic_int);
    FOR (idir)
    {
        if (is_periodic_int[idir] != 0 && is_periodic_int[idir] != 1)
        {
            geom_pp.error("is_periodic", "entries must be either 0 or 1");
        }
    }

    const auto check_boundary_conditions =
        [&boundary_pp,
         &is_periodic_int](const char *a_name,
                           const std::array<int, AMREX_SPACEDIM> &a_conditions)
    {
        std::string ignored_directions;
        std::string unset_directions;
        FOR (idir)
        {
            if (is_periodic_int[idir] == 1 && a_conditions[idir] != UNSET_BC)
            {
                if (!ignored_directions.empty())
                {
                    ignored_directions += ", ";
                }
                ignored_directions += static_cast<char>('x' + idir);
            }
            else if (is_periodic_int[idir] == 0 &&
                     a_conditions[idir] == UNSET_BC)
            {
                if (!unset_directions.empty())
                {
                    unset_directions += ", ";
                }
                unset_directions += static_cast<char>('x' + idir);
            }
        }

        if (!unset_directions.empty())
        {
            boundary_pp.error(
                a_name,
                "entries must be specified for the following non-periodic "
                "directions: " +
                    unset_directions);
        }

        if (!ignored_directions.empty())
        {
            boundary_pp.warning(
                a_name,
                "entries for the following periodic directions are ignored: " +
                    ignored_directions);
        }
    };

    check_boundary_conditions("hi_condition", hi_conditions);
    check_boundary_conditions("lo_condition", lo_conditions);
}

/// Set the geometry-dependent members and load the boundary parameters
void BoundaryConditions::define(const amrex::Geometry &a_geom)
{
    m_params.fill_params();
    GRParmParse pp;
    pp.get("evolution.num_ghosts", m_num_ghosts);

    pp.getarr("geometry.center", m_center);
    m_geom = a_geom;
}

/// change the asymptotic values of the variables for the Sommerfeld BCs
/// this will allow them to evolve during a simulation if necessary
void BoundaryConditions::set_vars_asymptotic_values(
    std::array<amrex::Real, NUM_VARS> &vars_asymptotic_values)
{
    m_params.vars_asymptotic_values = vars_asymptotic_values;
    m_asymptotic_values.clear();
}

/// The function which returns the parity of each of the vars in
/// StateVariables.hpp (It is only required for reflective boundary conditions.)
// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
int BoundaryConditions::get_state_var_parity(int a_comp, int a_dir)
{
    BCParity comp_parity = StateVariables::parities[a_comp];

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        comp_parity != BCParity::undefined,
        "BoundaryConditions: cannot have undefined variable parities if using "
        "reflective BCs");

    auto dir_parities = bc_parity_map.at(comp_parity);

    return dir_parities[a_dir];
}

/// Get the boundary condition for given face
int BoundaryConditions::get_boundary_condition(amrex::Orientation face) const
{
    return face.isLow() ? m_params.lo_condition[face.coordDir()]
                        : m_params.hi_condition[face.coordDir()];
}

bool BoundaryConditions::is_periodic(const int a_dir) const
{
    return m_params.is_periodic[a_dir];
}

void BoundaryConditions::apply_sommerfeld_boundaries(
    amrex::MultiFab &a_rhs, const amrex::MultiFab &a_soln) const
{
    // First-order extrapolation and reflection are imposed by AMReX when it
    // fills physical ghost cells. Sommerfeld instead requires this additional
    // replacement of the evolution RHS in the outer valid cells.
    if (!m_params.boundary_exists(SOMMERFELD_BC))
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
        {
            // These faces have usable ghost cells, so the central derivative
            // stencil can extend beyond the valid domain.
            domain.grow(face);
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
                                d1 = (0.5 / dx) *
                                     (sol(iv_offset1, n) - sol(iv_offset2, n));
                            }
                            // for each direction add dphidx * x^i
                            tmp += -d1 * loc[idir2];
                        }
                        // asymptotic values - these need to have been set in
                        // the params file
                        amrex::Real radius =
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
