#include "KleinGordonLevel.hpp"
#include <numeric>

using namespace amrex;

amrex::Vector<std::string> KleinGordonLevel::diagnostics;

void KleinGordonLevel::variableSetUp()
{
    BL_PROFILE("KleinGordonLevel::variableSetUp()");

    // Set up the state variables
    stateVariableSetUp();

    // Set up derived variables
    derive_lst.add(
        "analytic_soln", amrex::IndexType::TheCellType(), 1, diagnostics,
        derive_func_fab,
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

    constexpr Real t0 = -5.4;

    amrex::Vector<amrex::Real> start_times{-5.4, 5.4};
    amrex::Vector<amrex::Real> start_pos{
        midpts[0], midpts[1], midpts[2] + 0.5 * midpts[2],
        midpts[0], midpts[1], midpts[2] - 0.5 * midpts[2]};

    InitialConditions SineGordon(simParams().alpha, simParams().k_r);

    amrex::ParallelFor(
        S_new,
        [=] AMREX_GPU_DEVICE(int bi, int i, int j, int k) noexcept
        {
            Real x = problo[0] + (i + 0.5) * dx[0];
            Real y = problo[1] + (j + 0.5) * dx[1];
            Real z = problo[2] + (k + 0.5) * dx[2];

            snew[bi](i, j, k, 0) =
                SineGordon.breather_solution(x - midpts[0], 0);
            snew[bi](i, j, k, 1) =
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

void KleinGordonLevel::specificEvalRHS(amrex::MultiFab &a_soln,
                                       amrex::MultiFab &a_rhs,
                                       const double a_time)
{
    BL_PROFILE("KleinGordonLevel::specificEvalRHS()");

    const auto dxinv = Geom().InvCellSizeArray();
    AMREX_D_TERM(Real dx2inv = dxinv[0] * dxinv[0];
                 , Real dy2inv = dxinv[1] * dxinv[1];
                 , Real dz2inv = dxinv[2] * dxinv[2]);
    auto const &sa   = a_soln.arrays();
    auto const &sdot = a_rhs.arrays();

    Potential my_potential(simParams().scalar_mass);

    amrex::ParallelFor(
        a_soln,
        [=] AMREX_GPU_DEVICE(int bi, int i, int j, int k) noexcept
        {
            auto const &s = sa[bi];
            auto const &f = sdot[bi];

            //      Real phi2 = std::pow(s(i,j,k,0),2)+std::pow(s(i,j,k,2),2);
            amrex::Real phi = 0;

            f(i, j, k, 0) = s(i, j, k, 1);

            AMREX_D_TERM(
                Real lapx =
                    dx2inv *
                    (-2.5 * s(i, j, k, 0) +
                     (4. / 3.) * (s(i - 1, j, k, 0) + s(i + 1, j, k, 0)) -
                     (1. / 12.) * (s(i - 2, j, k, 0) + s(i + 2, j, k, 0)));
                , Real lapy =
                      dy2inv *
                      (-2.5 * s(i, j, k, 0) +
                       (4. / 3.) * (s(i, j - 1, k, 0) + s(i, j + 1, k, 0)) -
                       (1. / 12.) * (s(i, j - 2, k, 0) + s(i, j + 2, k, 0)));
                , Real lapz =
                      dz2inv *
                      (-2.5 * s(i, j, k, 0) +
                       (4. / 3.) * (s(i, j, k - 1, 0) + s(i, j, k + 1, 0)) -
                       (1. / 12.) * (s(i, j, k - 2, 0) + s(i, j, k + 2, 0))));

            f(i, j, k, 1) = AMREX_D_TERM(lapx, +lapy, +lapz);

            f(i, j, k, 1) -= std::sin(s(i, j, k, 0));
        });
    Gpu::streamSynchronize();
}

void KleinGordonLevel::errorEst(TagBoxArray &tags, int /*clearval*/,
                                int /*tagval*/, Real /*time*/,
                                int /*n_error_buf*/, int /*ngrow*/)
{
    BL_PROFILE("KleinGordonLevel::errorEst()");

    auto const &S_new = get_new_data(State_Type);

    const char tagval = TagBox::SET;
    auto const &a     = tags.arrays();
    auto const &s     = S_new.const_arrays();
    amrex::ParallelFor(tags,
                       [=] AMREX_GPU_DEVICE(int bi, int i, int j, int k)
                       {
                           // Just an example, not necessarily good choice.
                           if (amrex::Math::abs(s[bi](i, j, k, 1)) > 0.75)
                           {
                               a[bi](i, j, k) = tagval;
                           }
                       });
}
