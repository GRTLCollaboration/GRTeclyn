#include "KleinGordonLevel.hpp"
#include "FixedGridsTagger.hpp"
#include "FourthOrderDerivatives.hpp"
#include "KleinGordonRHS.hpp"
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
    derive_lst.addComponent("analytic_soln", desc_lst, State_Type, 0, 1);

    // The following is an example of how to use the current state to compute a
    // new derived variable that depends on the state variables and the
    // potential

    amrex::ParmParse pp("klein_gordon");
    std::string model{};
    pp.query("model", model);

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

    derive_lst.addComponent("rho", desc_lst, State_Type, 0, NUM_VARS);
}

void KleinGordonLevel::initData()
{
    BL_PROFILE("KleinGordonLevel::initData()");

    std::array<double, AMREX_SPACEDIM> center{};
    std::string model{};
    amrex::Real initial_time{0.0};

    amrex::ParmParse pp;
    pp.query("center", center);
    pp.query("klein_gordon.model", model);
    pp.query("klein_gordon.initial_time", initial_time);

    amrex::MultiFab &state_new = get_new_data(State_Type);
    auto const &array_new      = state_new.arrays();

    int dcomp{0};
    const amrex::Real current_time{
        0.0}; // initial time is an internal parameter
              // to the model class so the actual
              // simulation time is what we want here

    // NB: the analytic solutions are defined in InitialConditions.cpp
    // The functions below are defined in DerivedVariables.cpp
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
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
void KleinGordonLevel::specificEvalRHS(amrex::MultiFab &a_soln,
                                       amrex::MultiFab &a_rhs,
                                       const double a_time)
{
    BL_PROFILE("KleinGordonLevel::specificEvalRHS()");

    const auto dx         = Geom().CellSize(0);
    auto const &soln_arrs = a_soln.const_arrays();
    auto const &rhs_arrs  = a_rhs.arrays();

    amrex::ParmParse pp("klein_gordon");

    std::string model{};

    pp.query("model", model);

    // Here std::variant is used to allow dynamic polymorphism between two
    // unrelated classes. my_model_variant can be either Wave or SineGordon,
    // until the if statement selects a class based on the value of the string,
    // model.

    using ModelVariant = std::variant<Wave, SineGordon>;
    ModelVariant my_model_variant;

    if (model == "Wave")
    {
        my_model_variant = Wave{};
    }
    else
    {
        my_model_variant = SineGordon{};
    }

    // But! The KleinGordonRHS constructor doesn't accept a std::variant object
    // in it's constructor. So std::visit is called to handle the multiple types
    // that might be contained in a std::variant object. In this case, they are
    // both model_t type but a variant could also be constructed from
    // std::variant<string, int> or other such combo.

    // Points to notice:
    // * std::visit is using another lambda function as a visitor function on
    // top of the one in ParallelFor (but captured by reference this time [&] )
    // * the double ampersand is an rvalue reference - it prevents another copy
    // being made since the object can be initialized by the move constructor
    // instead

    std::visit(
        [&](auto &&my_model)
        {
            KleinGordonRHS rhs(simParams().sigma, dx, my_model);
            amrex::ParallelFor(
                a_soln,
                [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
                { rhs.compute(i, j, k, soln_arrs[box_no], rhs_arrs[box_no]); });
        },
        my_model_variant);

    amrex::Gpu::streamSynchronize();
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

    amrex::MultiFab &state_new = get_new_data(State_Type);

    const auto &tag_arrs   = tags.arrays();
    const auto &state_arrs = state_new.arrays();

    const amrex::Real dx         = Geom().CellSize(0);
    const int current_level      = Level();
    const amrex::Real box_length = Geom().ProbLength(0);
    std::array<double, AMREX_SPACEDIM> center{AMREX_D_DECL(0., 0., 0.)};
    GRParmParse pp;
    pp.query("center", center);

    FixedGridsTagger my_tagging_criterion{dx, current_level, box_length,
                                          center};

    amrex::ParallelFor(tags,
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       { my_tagging_criterion(i, j, k, tag_arrs[box_no]); });
    amrex::Gpu::streamSynchronize();
}
