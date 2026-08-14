/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#include "KleinGordonLevel.hpp"
#include "FixedGridsTagger.hpp"
#include "FourthOrderDerivatives.hpp"
#include "KleinGordonRHS.hpp"
#include "StateTypes.hpp"
#include <numeric>

void KleinGordonLevel::variableSetUp()
{
    BL_PROFILE("KleinGordonLevel::variableSetUp()");

    // Set up the state variables
    stateVariableSetUp();

    // The first two derived variables calculate the analytic solution
    //  for phi and Pi

    const std::string &comp_type                = {"analytic_soln"};
    const amrex::Vector<std::string> comp_names = {"phi_analytic",
                                                   "Pi_analytic"};

    int ncomp_analytic{
        static_cast<int>(comp_names.size())}; // how many derived variables

    derive_lst.add(
        comp_type, amrex::IndexType::TheCellType(), ncomp_analytic, comp_names,
        calc_analytic_solution, [=](const amrex::Box &box) { return box; },
        &amrex::cell_quartic_interp);
    derive_lst.addComponent("analytic_soln", desc_lst, state_index, 0, 1);

    // The following is an example of how to use the current state to compute a
    // new derived variable that depends on the state variables and the
    // potential

    amrex::ParmParse pp("klein_gordon");
    std::string model{};
    pp.get("model", model);

    const int ncomp_rho{1}; // only one component associated with energy density
    const int nghosts_rho{2};

    if (model == "Wave")
    {
        derive_lst.add(
            "rho", amrex::IndexType::TheCellType(), ncomp_rho,
            calc_energy_density<Wave>, [=](const amrex::Box &box)
            { return amrex::grow(box, nghosts_rho); },
            &amrex::cell_quartic_interp);
    }

    if (model.find("SineGordon") == 0)
    {
        derive_lst.add(
            "rho", amrex::IndexType::TheCellType(), ncomp_rho,
            calc_energy_density<SineGordon>, [=](const amrex::Box &box)
            { return amrex::grow(box, nghosts_rho); },
            &amrex::cell_quartic_interp);
    }

    derive_lst.addComponent("rho", desc_lst, state_index, 0, NUM_VARS);
}

void KleinGordonLevel::initData()
{
    BL_PROFILE("KleinGordonLevel::initData()");

    std::array<double, AMREX_SPACEDIM> center{};
    std::string model{};
    amrex::Real initial_time{0.0};

    amrex::ParmParse pp;
    pp.get("grteclyn.center", center);
    pp.get("klein_gordon.model", model);
    pp.get("klein_gordon.initial_time", initial_time);

    amrex::MultiFab &state_new = get_new_data(state_index);
    auto const &array_new      = state_new.arrays();

    int dcomp{0};
    const amrex::Real current_time{
        0.0}; // initial time is an internal parameter
              // to the model class so the actual
              // simulation time is what we want here

    // NB: the analytic solutions are defined in InitialConditions.cpp
    // The functions below are defined in DerivedVariables.cpp
    // NOLINTBEGIN(bugprone-branch-clone)
    if (model == "Wave")
    {
        calc_analytic_mf_3d<Wave>(state_new, dcomp, geom, current_time);
    }
    else if (model == "SineGordon1D")
    {
        calc_analytic_mf_1d<SineGordon>(state_new, dcomp, geom, current_time);
    }
    else
    {
        calc_analytic_mf_3d<SineGordon>(state_new, dcomp, geom, current_time);
    }
    // NOLINTEND(bugprone-branch-clone)
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
void KleinGordonLevel::specificEvalRHS(amrex::MultiFab &a_soln,
                                       amrex::MultiFab &a_rhs,
                                       const double a_time)
{
    BL_PROFILE("KleinGordonLevel::specificEvalRHS()");

    amrex::ParmParse pp("klein_gordon");

    std::string model{};

    pp.query("model", model);

    if (model == "Wave")
    {
        eval_model_specific_rhs<Wave>(a_soln, a_rhs);
    }
    else
    {
        eval_model_specific_rhs<SineGordon>(a_soln, a_rhs);
    }

    amrex::Gpu::streamSynchronize();
}

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
template <class model_t>
void KleinGordonLevel::eval_model_specific_rhs(amrex::MultiFab &a_soln,
                                               amrex::MultiFab &a_rhs)
// NOLINTEND(bugprone-easily-swappable-parameters)
{

    const auto dx                 = Geom().CellSize(0);
    const auto &const_soln_arrays = a_soln.const_arrays();
    const auto &rhs_arrays        = a_rhs.arrays();

    model_t my_model;
    KleinGordonRHS kg_rhs(dx, my_model);

    amrex::ParallelFor(
        a_soln,
        [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz) noexcept
        { kg_rhs(ix, iy, iz, rhs_arrays[box_no], const_soln_arrays[box_no]); });
}

void KleinGordonLevel::tag_cells(amrex::TagBoxArray &tags,
                                 amrex::Real a_regrid_threshold)
{

    // This is an example of how to apply a tagging criterion to mark cells for
    // further refinement.

    // The FixedGridsTagger used here will always tag the inner
    // L/4 cells of the simulation volume, regardless of the
    // field values. You can also ask that the cells be tagged
    // when the field value or its derivative reaches a certain value.

    BL_PROFILE("KleinGordonLevel::tag_cells()");

    amrex::MultiFab &state_new = get_new_data(state_index);

    const auto &tag_arrs   = tags.arrays();
    const auto &state_arrs = state_new.arrays();

    const amrex::Real dx         = Geom().CellSize(0);
    const int current_level      = Level();
    const amrex::Real box_length = Geom().ProbLength(0);
    std::array<double, AMREX_SPACEDIM> center{AMREX_D_DECL(0., 0., 0.)};
    GRParmParse pp;
    pp.query("grteclyn.center", center);

    FixedGridsTagger my_tagging_criterion{dx, current_level, box_length,
                                          center};

    amrex::ParallelFor(tags,
                       [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
                       { my_tagging_criterion(ix, iy, iz, tag_arrs[box_no]); });
    amrex::Gpu::streamSynchronize();
}
