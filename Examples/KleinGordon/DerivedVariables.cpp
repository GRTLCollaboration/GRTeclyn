#include "DerivedVariables.hpp"

void calc_analytic_solution_mf(amrex::MultiFab &mf_out, int dcomp,
                               int /*numcomp*/, const amrex::MultiFab &mf_in,
                               const amrex::Geometry &geom,
                               const amrex::Real time, const int * /*bcomp*/,
                               int /*scomp*/)

{
    amrex::ParmParse pp;

    // Which model is it?
    std::string model;
    pp.query("model", model);

    // Get the geometry of the problem
    const auto problo = geom.ProbLoArray();
    const auto dx     = geom.CellSizeArray();

    std::array<double, AMREX_SPACEDIM> center{};
    pp.query("center", center);

    auto const &array4_out = mf_out.arrays();

    if (model == "SineGordon1D")
    {
        amrex::Real alpha = 0.7;
        pp.query("alpha", alpha);

        InitialConditions AnalyticSoln(alpha);

        amrex::ParallelFor(
            mf_out, mf_out.nGrowVect(),
            [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
            {
                amrex::Real x = problo[0] + (i + 0.5) * dx[0] - center[0];

                amrex::Real exact_soln =
                    AnalyticSoln.breather_solution(x, time);

                array4_out[box_no](i, j, k, dcomp) = exact_soln;
            });
    }

    if (model == "SineGordon3D")
    {
        amrex::Real alpha        = 0.7;
        amrex::Real initial_time = -5.4;

        pp.query("alpha", alpha);
        pp.query("initial_time", initial_time);

        InitialConditions AnalyticSoln(alpha);

        amrex::ParallelFor(
            mf_out, mf_out.nGrowVect(),
            [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
            {
                amrex::Real x = problo[0] + (i + 0.5) * dx[0] - center[0];
                amrex::Real y = problo[1] + (j + 0.5) * dx[1] - center[1];
                amrex::Real z = problo[2] + (k + 0.5) * dx[2] - center[2];

                array4_out[box_no](i, j, k, dcomp) =
                    AnalyticSoln.breather_solution(x, y, z,
                                                   initial_time + time);
            });
    }

    if (model == "Wave")
    {
        amrex::Real k_r = 1;
        pp.query("wave_vector", k_r);

        InitialConditions AnalyticSoln(k_r);

        amrex::ParallelFor(
            mf_out, mf_out.nGrowVect(),
            [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
            {
                amrex::Real x = problo[0] + (i + 0.5) * dx[0] - center[0];
                amrex::Real y = problo[1] + (j + 0.5) * dx[1] - center[1];
                amrex::Real z = problo[2] + (k + 0.5) * dx[2] - center[2];

                amrex::Real exact_soln =
                    AnalyticSoln.travelling_wave(x, y, z, time);

                array4_out[box_no](i, j, k, dcomp) = exact_soln;
            });
    }
}
