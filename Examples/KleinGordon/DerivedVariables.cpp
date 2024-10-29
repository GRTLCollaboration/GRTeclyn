#include <DerivedVariables.hpp>

void derive_func_fab(const amrex::Box &bx, amrex::FArrayBox &derfab, int dcomp,
                     int /*numcomp*/, const amrex::FArrayBox &datfab,
                     const amrex::Geometry &geom, const amrex::Real time,
                     const int * /*bcomp*/, int /*scomp*/)

{
    amrex::ParmParse pp;

    amrex::Real k_r = 1;
    pp.query("wave_vector", k_r);

    amrex::Real alpha = 0.7;
    pp.query("alpha", alpha);

    InitialConditions SineGordon(alpha, k_r);

    const auto problo = geom.ProbLoArray();
    const auto probhi = geom.ProbHiArray();
    const auto dx     = geom.CellSizeArray();

    amrex::Real midpts[3];
    midpts[0] = 0.5 * (probhi[0] - problo[0]);
    midpts[1] = 0.5 * (probhi[1] - problo[1]);
    midpts[2] = 0.5 * (probhi[2] - problo[2]);

    auto const &s     = datfab.array();
    auto const &s_out = derfab.array();

    amrex::ParallelFor(bx,
                       [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                       {
                           amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                           amrex::Real y = problo[1] + (j + 0.5) * dx[1];
                           amrex::Real z = problo[2] + (k + 0.5) * dx[2];

                           amrex::Real exact_soln =
                               SineGordon.breather_solution(x - midpts[0],
                                                            time);

                           s_out(i, j, k, dcomp) = exact_soln;
                       });
}
