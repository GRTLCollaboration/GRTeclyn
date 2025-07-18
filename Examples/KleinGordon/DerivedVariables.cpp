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
        amrex::Real alpha{0.7};
        pp.query("alpha", alpha);
        SineGordon my_sine_gordon_model(alpha, initial_time);
        my_sine_gordon_model.calc_mf_1d(mf_out, dcomp, geom, time);
    }
    if (model == "SineGordon3D")
    {
        amrex::Real alpha{0.7};
        pp.query("alpha", alpha);
        SineGordon my_sine_gordon_model(alpha, initial_time);
        my_sine_gordon_model.calc_mf_3d(mf_out, dcomp, geom, time);
    }
    // This is a special case because the analytic solution assumes no potential
    // If a potential is given then fill the analytic solution with zeros
    if (model == "Wave" && scalar_mass == 0)
    {
        amrex::Real k_r{1.0};
        pp.query("wave_vector", k_r);
        Wave my_wave_model(k_r, scalar_mass, initial_time);
        my_wave_model.calc_mf(mf_out, dcomp, geom, time);
    }
}
