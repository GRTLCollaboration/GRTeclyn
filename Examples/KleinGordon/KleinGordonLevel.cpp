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

    // Set up derived variables
    derive_lst.add(
        "analytic_soln", amrex::IndexType::TheCellType(), 1, calc_derive_mf,
        [=](const amrex::Box &box)
        { return amrex::grow(box, simParams().num_ghosts); },
        &amrex::cell_quartic_interp);

    derive_lst.addComponent("analytic_soln", desc_lst, State_Type, 0, 1);
}

void KleinGordonLevel::initData()
{
    BL_PROFILE("KleinGordonLevel::initData()");

    const auto problo = geom.ProbLoArray();
    const auto dx     = geom.CellSizeArray();

    std::array<double, AMREX_SPACEDIM> center{};

    amrex::ParmParse pp;
    pp.query("center", center);

    MultiFab &S_new  = get_new_data(State_Type);
    auto const &snew = S_new.arrays();

    if (simParams().model == "Wave")
    {
        InitialConditions Wave(simParams().k_r);
        amrex::ParallelFor(
            S_new,
            [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
            {
                amrex::Real x = problo[0] + (i + 0.5) * dx[0] - center[0];
                amrex::Real y = problo[1] + (j + 0.5) * dx[1] - center[1];
                amrex::Real z = problo[2] + (k + 0.5) * dx[2] - center[2];

                snew[box_no](i, j, k, 0) = Wave.travelling_wave(x, y, z, 0);
                snew[box_no](i, j, k, 1) =
                    Wave.travelling_wave_deriv(x, y, z, 0);
            });
    }
    else
    {
        InitialConditions SineGordon(simParams().alpha);

        if (simParams().model == "SineGordon1D")
        {
            amrex::ParallelFor(
                S_new,
                [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
                {

                    amrex::Real x = problo[0] + (i + 0.5) * dx[0] - center[0];

                    snew[box_no](i, j, k, 0) =
                        SineGordon.breather_solution(x, 0);
                    snew[box_no](i, j, k, 1) =
                        SineGordon.breather_solution_deriv(x, 0);
                });
        }
        else
        {

            amrex::Real initial_time = -5.4;

            amrex::ParmParse pp;
            pp.query("initial_time", initial_time);

            amrex::ParallelFor(
                S_new,
                [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
                {

                    amrex::Real x = problo[0] + (i + 0.5) * dx[0] - center[0];
                    amrex::Real y = problo[1] + (j + 0.5) * dx[1] - center[1];
                    amrex::Real z = problo[2] + (k + 0.5) * dx[2] - center[2];

                    snew[box_no](i, j, k, 0) =
                        SineGordon.breather_solution(x, y, z, initial_time);
                    snew[box_no](i, j, k, 1) = 0;
                });
        }
    }
    amrex::Gpu::streamSynchronize();
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
