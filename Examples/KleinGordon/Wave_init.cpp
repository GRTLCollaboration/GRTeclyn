#include "InitialConditions.H"
#include "KleinGordonLevel.hpp"
#include <cmath>

using namespace amrex;

void KleinGordonLevel::initData()
{
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

    InitialConditions SineGordon(alpha, k_r);

    amrex::ParallelFor(
        S_new,
        [=] AMREX_GPU_DEVICE(int bi, int i, int j, int k) noexcept
        {
            Real x = problo[0] + (i + 0.5) * dx[0];
            Real y = problo[1] + (j + 0.5) * dx[1];
            Real z = problo[2] + (k + 0.5) * dx[2];

            //	Real rr2 = (x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)+(z-0.5)*(z-0.5);

            //        constexpr Real Pi = 3.1415926535897932384626;
            // constexpr Real mu = 0.7;
            // constexpr Real mu_coeff = mu/(std::sqrt(1-mu*mu));
            // constexpr Real t0 = 0;

            // constexpr Real k_r = 1;
            // constexpr Real omega = 1;
            for (int n = 0; n < nfields; n++)
            {
                //	    snew[bi](i,j,k,2*n) = 0.0;
                //	    snew[bi](i,j,k,2*n+1) = std::exp(-16.*rr2) *
                // std::pow(std::cos(Pi*rr2),6);

                // snew[bi](i,j,k,2*n) = 1+ampl[n]*std::exp(-width[n]*rr2);
                // snew[bi](i,j,k,2*n+1) = 0;

                snew[bi](i, j, k, 2 * n) =
                    SineGordon.breather_solution(x - midpts[0], 0);
                snew[bi](i, j, k, 2 * n + 1) =
                    SineGordon.breather_solution_deriv(x - midpts[0], 0);

                //	    snew[bi](i,j,k,2*n) =
                // SineGordon.breather_solution(x-start_pos[0], y-start_pos[1],
                // z-start_pos[2], start_times[0]) +
                // SineGordon.breather_solution(x-start_pos[3], y-start_pos[4],
                // z-start_pos[5], start_times[1]); snew[bi](i,j,k,2*n+1) = 0;
            }
        });
}
