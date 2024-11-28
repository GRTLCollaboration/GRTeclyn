#include "DerivedVariables.hpp"

void calc_derive_mf(amrex::MultiFab &mf_out, int dcomp, int /*numcomp*/,
                    const amrex::MultiFab &mf_in, const amrex::Geometry &geom,
                    const amrex::Real time, const int * /*bcomp*/,
                    int /*scomp*/)

{
    amrex::ParmParse pp;

    // Which model is it?
    std::string model;
    pp.query("model", model);

    // Get the geometry of the problem
    const auto problo = geom.ProbLoArray();
    const auto probhi = geom.ProbHiArray();
    const auto dx     = geom.CellSizeArray();

    amrex::Real midpts[3];

    int i = 0;

    FOR (i)
        midpts[i] = 0.5 * (probhi[i] - problo[i]);

    auto const &array4_out = mf_out.arrays();

    if (model == "SineGordon")
    {
        amrex::Real alpha = 0.7;
        pp.query("alpha", alpha);

        InitialConditions AnalyticSoln(alpha);

        amrex::ParallelFor(
            mf_out, mf_out.nGrowVect(),
            [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
            {
                amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                amrex::Real y = problo[1] + (j + 0.5) * dx[1];
                amrex::Real z = problo[2] + (k + 0.5) * dx[2];

                amrex::Real exact_soln =
                    AnalyticSoln.breather_solution(x - midpts[0], time);

                array4_out[box_no](i, j, k, dcomp) = exact_soln;
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
                amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                amrex::Real y = problo[1] + (j + 0.5) * dx[1];
                amrex::Real z = problo[2] + (k + 0.5) * dx[2];

                amrex::Real exact_soln = AnalyticSoln.travelling_wave(
                    x - midpts[0], y - midpts[1], z - midpts[2], time);

                array4_out[box_no](i, j, k, dcomp) = exact_soln;
            });
    }
}
