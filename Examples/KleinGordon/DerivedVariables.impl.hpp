#if !defined(DERIVEDVARIABLES_HPP_)
#error "This file should only be included through DerivedVariables.hpp"
#endif

#ifndef DERIVEDVARIABLES_IMPL_HPP_
#define DERIVEDVARIABLES_IMPL_HPP_

void AMREX_FORCE_INLINE calc_analytic_solution(
    amrex::MultiFab &mf_out, int dcomp, int /*numcomp*/,
    const amrex::MultiFab & /*mf_in*/, const amrex::Geometry &geom,
    const amrex::Real time, const int * /*bcomp*/, int /*scomp*/)
{
    amrex::ParmParse pp;
    std::string model{};
    pp.query("model", model);

    amrex::Real scalar_mass{0.0};
    pp.query("scalar_mass", scalar_mass);

    if (model == "SineGordon1D")
    {
        calc_analytic_mf_1d<SineGordon>(mf_out, dcomp, geom, time);
    }
    if (model == "SineGordon3D")
    {
        calc_analytic_mf_3d<SineGordon>(mf_out, dcomp, geom, time);
    }
    // This is a special case because the analytic solution assumes no potential
    // If a potential is given then fill the analytic solution with zeros
    if (model == "Wave" && scalar_mass == 0)
    {
        calc_analytic_mf_3d<Wave>(mf_out, dcomp, geom, time);
    }
}

template <typename model_t>
AMREX_FORCE_INLINE void calc_analytic_mf_3d(amrex::MultiFab &mf_out, int dcomp,
                                            const amrex::Geometry &geom,
                                            const amrex::Real time)
{
    amrex::ParmParse pp;

    // Get the geometry of the simulation
    const auto dx = geom.CellSizeArray();

    std::array<amrex::Real, AMREX_SPACEDIM> center{};
    pp.query("center", center);

    model_t model;

    auto const &arrs = mf_out.arrays();

    amrex::ParallelFor(
        mf_out, mf_out.nGrowVect(),
        [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
        {
            amrex::IntVect grid_pos(i, j, k);
            Coordinates pos(grid_pos, dx[0], center);

            arrs[box_no](i, j, k, dcomp) =
                model.calculate(pos.x, pos.y, pos.z, time);
            arrs[box_no](i, j, k, dcomp + 1) =
                model.calculate_derivative(x, y, z, time);
        });
    amrex::Gpu::streamSynchronize();
};

template <typename model_t>
AMREX_FORCE_INLINE void calc_analytic_mf_1d(amrex::MultiFab &mf_out, int dcomp,
                                            const amrex::Geometry &geom,
                                            const amrex::Real time)
{
    amrex::ParmParse pp;

    // Get the geometry of the simulation
    const auto dx = geom.CellSizeArray();

    std::array<amrex::Real, AMREX_SPACEDIM> center{};
    pp.query("center", center);

    model_t model;

    auto const &arrs = mf_out.arrays();

    amrex::ParallelFor(
        mf_out, mf_out.nGrowVect(),
        [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
        {
            amrex::IntVect grid_pos(i, j, k);
            Coordinates pos(grid_pos, dx[0], center);

            arrs[box_no](i, j, k, dcomp) = model.calculate(pos.x, time);
            arrs[box_no](i, j, k, dcomp + 1) =
                model.calculate_derivative(x, time);
        });
    amrex::Gpu::streamSynchronize();
};

#endif // DERIVEDVARIABLES_IMPL_HPP_
