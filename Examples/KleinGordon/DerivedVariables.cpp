#include "DerivedVariables.hpp"

void calc_sine_gordon_1d_analytic_solution(amrex::MultiFab &mf_out, int dcomp,
                                           int /*numcomp*/,
                                           const amrex::MultiFab &mf_in,
                                           const amrex::Geometry &geom,
                                           const amrex::Real time,
                                           const int * /*bcomp*/, int /*scomp*/)
{

    // Get the geometry of the problem
    const auto problo = geom.ProbLoArray();
    const auto dx     = geom.CellSizeArray();

    amrex::ParmParse pp;

    std::array<double, AMREX_SPACEDIM> center{};
    pp.query("center", center);

    amrex::Real alpha = 0.7;
    pp.query("alpha", alpha);

    auto const &arrs = mf_out.arrays();

    amrex::ParallelFor(
        mf_out, mf_out.nGrowVect(),
        [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
        {
            amrex::Real x = problo[0] + (i + 0.5) * dx[0] - center[0];

            amrex::Real exact_soln =
                KleinGordon::breather_solution(alpha, x, time);

            arrs[box_no](i, j, k, dcomp) = exact_soln;
        });
}

void calc_sine_gordon_3d_analytic_solution(amrex::MultiFab &mf_out, int dcomp,
                                           int /*numcomp*/,
                                           const amrex::MultiFab &mf_in,
                                           const amrex::Geometry &geom,
                                           const amrex::Real time,
                                           const int * /*bcomp*/, int /*scomp*/)

{

    // Get the geometry of the problem
    const auto problo = geom.ProbLoArray();
    const auto dx     = geom.CellSizeArray();

    amrex::ParmParse pp;

    std::array<double, AMREX_SPACEDIM> center{};
    pp.query("center", center);

    amrex::Real alpha = 0.7;
    pp.query("alpha", alpha);

    amrex::Real initial_time = -5.4;
    pp.query("initial_time", initial_time);

    auto const &arrs = mf_out.arrays();

    amrex::ParallelFor(
        mf_out, mf_out.nGrowVect(),
        [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
        {
            amrex::Real x = problo[0] + (i + 0.5) * dx[0] - center[0];
            amrex::Real y = problo[1] + (j + 0.5) * dx[1] - center[1];
            amrex::Real z = problo[2] + (k + 0.5) * dx[2] - center[2];

            arrs[box_no](i, j, k, dcomp) = KleinGordon::breather_solution(
                alpha, x, y, z, initial_time + time);
        });
}

void calc_wave_analytic_solution(amrex::MultiFab &mf_out, int dcomp,
                                 int /*numcomp*/, const amrex::MultiFab &mf_in,
                                 const amrex::Geometry &geom,
                                 const amrex::Real time, const int * /*bcomp*/,
                                 int /*scomp*/)

{
    amrex::ParmParse pp;

    // Get the geometry of the problem
    const auto problo = geom.ProbLoArray();
    const auto dx     = geom.CellSizeArray();

    std::array<double, AMREX_SPACEDIM> center{};
    pp.query("center", center);

    auto const &arrs = mf_out.arrays();

    amrex::Real k_r = 1;
    pp.query("wave_vector", k_r);

    amrex::ParallelFor(
        mf_out, mf_out.nGrowVect(),
        [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
        {
            amrex::Real x = problo[0] + (i + 0.5) * dx[0] - center[0];
            amrex::Real y = problo[1] + (j + 0.5) * dx[1] - center[1];
            amrex::Real z = problo[2] + (k + 0.5) * dx[2] - center[2];

            amrex::Real exact_soln =
                KleinGordon::travelling_wave(k_r, x, y, z, time);

            arrs[box_no](i, j, k, dcomp) = exact_soln;
        });
}
