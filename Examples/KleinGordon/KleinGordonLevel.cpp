#include "KleinGordonLevel.hpp"
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
    const auto probhi = geom.ProbHiArray();
    const auto dx     = geom.CellSizeArray();

    Real midpts[3];
    midpts[0] = 0.5 * (probhi[0] - problo[0]);
    midpts[1] = 0.5 * (probhi[1] - problo[1]);
    midpts[2] = 0.5 * (probhi[2] - problo[2]);

    MultiFab &S_new  = get_new_data(State_Type);
    auto const &snew = S_new.arrays();

    constexpr Real initial_time = -5.4;

    amrex::Vector<amrex::Real> start_times{initial_time, initial_time * -1.0};
    amrex::Vector<amrex::Real> start_pos{
        midpts[0], midpts[1], midpts[2] + 0.5 * midpts[2],
        midpts[0], midpts[1], midpts[2] - 0.5 * midpts[2]};

    InitialConditions SineGordon(simParams().alpha);

    amrex::ParallelFor(
        S_new,
        [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
        {
            Real x = problo[0] + (i + 0.5) * dx[0];
            Real y = problo[1] + (j + 0.5) * dx[1];
            Real z = problo[2] + (k + 0.5) * dx[2];

            snew[box_no](i, j, k, 0) =
                SineGordon.breather_solution(x - midpts[0], 0);
            snew[box_no](i, j, k, 1) =
                SineGordon.breather_solution_deriv(x - midpts[0], 0);

            //	    snew[bi](i,j,k,2*n) =
            // SineGordon.breather_solution(x-start_pos[0], y-start_pos[1],
            // z-start_pos[2], start_times[0]) +
            // SineGordon.breather_solution(x-start_pos[3], y-start_pos[4],
            // z-start_pos[5], start_times[1]); snew[bi](i,j,k,2*n+1) = 0;
        });
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

    const auto dx    = Geom().CellSize(0);
    const auto dxinv = Geom().InvCellSizeArray();
    AMREX_D_TERM(Real dx2inv = dxinv[0] * dxinv[0];
                 , Real dy2inv = dxinv[1] * dxinv[1];
                 , Real dz2inv = dxinv[2] * dxinv[2]);
    auto const &a_soln_a4 = a_soln.const_arrays();
    auto const &a_rhs_a4  = a_rhs.arrays();

    Potential my_potential(simParams().scalar_mass);

    KleinGordonRHS<FourthOrderDerivatives, Potential> klein_gordon_rhs(
        simParams().sigma, dx, my_potential);

    amrex::ParallelFor(
        a_soln,
        [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
            klein_gordon_rhs.compute(i, j, k, a_soln_a4[box_no],
                                     a_rhs_a4[box_no]);
        });

    Gpu::streamSynchronize();
}

void KleinGordonLevel::errorEst(TagBoxArray &tags, int /*clearval*/,
                                int /*tagval*/, Real /*time*/,
                                int /*n_error_buf*/, int /*ngrow*/)
{
    BL_PROFILE("KleinGordonLevel::errorEst()");

    auto const &state = get_new_data(State_Type);

    const char tagval     = TagBox::SET;
    auto const &arr       = tags.arrays();
    auto const &state_arr = state.const_arrays();
    amrex::ParallelFor(tags,
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       {
                           // Just an example, not necessarily good choice.
                           if (amrex::Math::abs(state_arr[box_no](i, j, k, 1)) >
                               0.75)
                           {
                               arr[box_no](i, j, k) = tagval;
                           }
                       });
}
