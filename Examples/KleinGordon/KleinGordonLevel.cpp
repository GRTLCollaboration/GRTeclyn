#include "KleinGordonLevel.hpp"
#include "FixedGridsTagger.hpp"
#include "FourthOrderDerivatives.hpp"
#include "KleinGordonRHS.hpp"
#include <numeric>

using namespace amrex;

void KleinGordonLevel::variableSetUp()
{
    BL_PROFILE("KleinGordonLevel::variableSetUp()");

    // Set up the state variables
    stateVariableSetUp();

    const std::string &comp_type                = {"analytic_soln"};
    const amrex::Vector<std::string> comp_names = {"phi_analytic",
                                                   "Pi_analytic"};

    derive_lst.add(
        comp_type, amrex::IndexType::TheCellType(),
        static_cast<int>(comp_names.size()), comp_names, calc_analytic_solution,
        [=](const amrex::Box &box)
        { return amrex::grow(box, simParams().num_ghosts); },
        &amrex::cell_quartic_interp);
    derive_lst.addComponent("analytic_soln", desc_lst, State_Type, 0, 1);
}

void KleinGordonLevel::initData()
{
    BL_PROFILE("KleinGordonLevel::initData()");

    std::array<double, AMREX_SPACEDIM> center{};
    std::string model{};

    amrex::ParmParse pp;
    pp.query("center", center);
    pp.query("model", model);

    MultiFab &state_new   = get_new_data(State_Type);
    auto const &array_new = state_new.arrays();

    int dcomp{0};
    amrex::Real initial_time{0.0};
    pp.query("initial_time", initial_time);

    // NB: the analytic solutions are defined in InitialConditions.cpp
    // The functions below are defined in DerivedVariables.cpp
    if (model == "Wave")
    {
        calc_wave_analytic_solution(state_new, dcomp, geom, initial_time);

        amrex::Real scalar_mass{0.0};
        pp.query("scalar_mass", scalar_mass);

        amrex::ParallelFor(
            state_new, state_new.nGrowVect(),
            [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
            {
                amrex::Real phi = array_new[box_no](i, j, k, dcomp);

                amrex::Real V_of_phi =
                    0.5 * scalar_mass * scalar_mass * phi * phi;

                amrex::Real dVdphi = scalar_mass * scalar_mass * phi;

                array_new[box_no](i, j, k, dcomp) += V_of_phi;

                array_new[box_no](i, j, k, dcomp + 1) += dVdphi;
            });
        amrex::Gpu::streamSynchronize();
    }
    else if (model == "SineGordon1D")
    {
        calc_sine_gordon_1d_analytic_solution(state_new, dcomp, geom,
                                              initial_time);
    }
    else
    {
        calc_sine_gordon_3d_analytic_solution(state_new, dcomp, geom,
                                              initial_time);
    }
}
void KleinGordonLevel::specificAdvance()
{
    // Usually constraints are enforced here, but we don't have any...
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

    Potential my_potential(simParams().scalar_mass);

    KleinGordonRHS<FourthOrderDerivatives, Potential> klein_gordon_rhs(
        simParams().sigma, dx, my_potential);

    amrex::ParallelFor(
        a_soln,
        [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
        {
            klein_gordon_rhs.compute(i, j, k, soln_arrs[box_no],
                                     rhs_arrs[box_no]);
        });

    amrex::Gpu::streamSynchronize();
}

void KleinGordonLevel::tag_cells(amrex::TagBoxArray &tags,
                                 amrex::Real a_regrid_threshold)
{
    BL_PROFILE("KleinGordonLevel::tag_cells()");

    amrex::MultiFab &state_new = get_new_data(State_Type);

    const char tagval      = TagBox::SET;
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
                       {
                           my_tagging_criterion.operator()<amrex::Real>(
                               i, j, k, tag_arrs[box_no]);
                       });
    amrex::Gpu::streamSynchronize();
}
