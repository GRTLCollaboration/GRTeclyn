#include "DerivedVariables.hpp"

void calc_analytic_solution(amrex::MultiFab &mf_out, int dcomp, int /*numcomp*/,
                            const amrex::MultiFab & /*mf_in*/,
                            const amrex::Geometry &geom, const amrex::Real time,
                            const int * /*bcomp*/, int /*scomp*/)
{
    amrex::ParmParse pp;
    std::string model{};
    pp.query("model", model);

    amrex::Real scalar_mass{0.0};
    pp.query("scalar_mass", scalar_mass);

    amrex::Real initial_time = -5.4;
    pp.query("initial_time", initial_time);

    if (model == "SineGordon1D")
    {
        calc_sine_gordon_1d_analytic_solution(mf_out, dcomp, geom,
                                              initial_time + time);
    }
    if (model == "SineGordon3D")
    {
        calc_sine_gordon_3d_analytic_solution(mf_out, dcomp, geom,
                                              initial_time + time);
    }
    // This is a special case because the analytic solution assumes no potential
    // If a potential is given then fill the analytic solution with zeros
    if (model == "Wave" && scalar_mass == 0)
    {
        calc_wave_analytic_solution(mf_out, dcomp, geom, initial_time + time);
    }
}

void calc_sine_gordon_1d_analytic_solution(amrex::MultiFab &mf_out, int dcomp,
                                           const amrex::Geometry &geom,
                                           const amrex::Real time)
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

            arrs[box_no](i, j, k, dcomp) =
                KleinGordon::breather_solution(alpha, x, time);

            arrs[box_no](i, j, k, dcomp + 1) =
                KleinGordon::breather_solution_deriv(alpha, x, time);
        });
    amrex::Gpu::streamSynchronize();
}

void calc_sine_gordon_3d_analytic_solution(amrex::MultiFab &mf_out, int dcomp,
                                           const amrex::Geometry &geom,
                                           const amrex::Real time)
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
            amrex::Real y = problo[1] + (j + 0.5) * dx[1] - center[1];
            amrex::Real z = problo[2] + (k + 0.5) * dx[2] - center[2];

            arrs[box_no](i, j, k, dcomp) =
                KleinGordon::breather_solution(alpha, x, y, z, time);
            arrs[box_no](i, j, k, dcomp + 1) = 0.0;
        });
    amrex::Gpu::streamSynchronize();
}

void calc_wave_analytic_solution(amrex::MultiFab &mf_out, int dcomp,
                                 const amrex::Geometry &geom,
                                 const amrex::Real time)
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

            arrs[box_no](i, j, k, dcomp) =
                KleinGordon::travelling_wave(k_r, x, y, z, time);
            arrs[box_no](i, j, k, dcomp + 1) =
                KleinGordon::travelling_wave_deriv(k_r, x, y, z, time);
        });
    amrex::Gpu::streamSynchronize();
}
