#ifndef RADIALRECIPE_MATTER_DISPATCH_HPP_
#define RADIALRECIPE_MATTER_DISPATCH_HPP_

#include "BiComplexScalarField.hpp"
#include "CCZ4RHSWithMatter.hpp"
#include "ComplexScalarField.hpp"
#include "ConstraintsWithMatter.hpp"
#include "ExoticScalarField.hpp"
#include "GRTresnaIndependentScalars.hpp"
#include "MovingPunctureGaugeWithMatter.hpp"
#include "RLMatterPumpParams.hpp"
#include "ScalarField.hpp"
#include "SimulationParameters.hpp"
#include "Weyl4WithMatter.hpp"

#include <AMReX_MultiFab.H>

namespace RadialRecipeMatter
{

inline bool uses_independent_scalars(const SimulationParameters &params)
{
    return params.recipe_matter_model == "grtresna_independent_scalars" &&
           params.recipe_num_scalar_fields > 0;
}

inline bool uses_complex_scalar(const SimulationParameters &params)
{
    return params.recipe_matter_model == "grtresna_complex_scalar";
}

inline bool uses_bicomplex_scalar(const SimulationParameters &params)
{
    return params.recipe_matter_model == "grtresna_bicomplex_scalar";
}

inline void setup_derived_quantities(int state_index,
                                     const SimulationParameters &params)
{
    if (uses_independent_scalars(params))
    {
        ConstraintsWithMatter<GRTresnaIndependentScalars>::set_up(state_index);
        Weyl4WithMatter<GRTresnaIndependentScalars>::set_up(state_index);
        return;
    }
    if (uses_complex_scalar(params))
    {
        ConstraintsWithMatter<ComplexScalarField>::set_up(state_index);
        Weyl4WithMatter<ComplexScalarField>::set_up(state_index);
        return;
    }
    if (uses_bicomplex_scalar(params))
    {
        ConstraintsWithMatter<BiComplexScalarField>::set_up(state_index);
        Weyl4WithMatter<BiComplexScalarField>::set_up(state_index);
        return;
    }
    if (params.recipe_exotic_matter)
    {
        ConstraintsWithMatter<ExoticScalarField<DefaultPotential>>::set_up(
            state_index);
        Weyl4WithMatter<ExoticScalarField<DefaultPotential>>::set_up(state_index);
        return;
    }
    ConstraintsWithMatter<ScalarField<DefaultPotential>>::set_up(state_index);
    Weyl4WithMatter<ScalarField<DefaultPotential>>::set_up(state_index);
}

//! Build the multi-site pump: spotlight centres follow the live tracked lump
//! centroids (RLRuntime::g_lump_state), per-lump amp/freq/phase come from the
//! (RL-updated) simulation parameters.  When RL is off every amplitude is 0, so
//! the pump is numerically inactive regardless of the (unseeded) centres.
inline RLMatterPumpParams build_rl_pump(const SimulationParameters &params,
                                        int num_fields)
{
    RLMatterPumpParams pump;
    int n = params.rl_num_lumps;
    if (n > GRTRESNA_MAX_INDEPENDENT_SCALARS)
        n = GRTRESNA_MAX_INDEPENDENT_SCALARS;
    if (n < 0)
        n = 0;
    pump.num_sites       = n;
    pump.width           = params.rl_pump_width;
    pump.governor_center = params.rl_l2_ham_governor_center;
    pump.governor_width  = params.rl_l2_ham_governor_width;
    pump.governor        = RLRuntime::tanh_governor(
        RLRuntime::cached_L2_Ham(), pump.governor_center, pump.governor_width);
    pump.num_fields      = num_fields;
    for (int s = 0; s < n; ++s)
    {
        pump.sites[s].center_x  = RLRuntime::g_lump_state[s].x;
        pump.sites[s].center_y  = RLRuntime::g_lump_state[s].y;
        pump.sites[s].center_z  = RLRuntime::g_lump_state[s].z;
        pump.sites[s].amplitude = params.rl_pump_amplitude[s];
        pump.sites[s].frequency = params.rl_pump_frequency[s];
        pump.sites[s].phase     = params.rl_pump_phase[s];
        // Route this spotlight to the canonical (+1) or phantom (-1) field for
        // the two-complex-field model.  Harmless for single-field models.
        pump.sites[s].field_sign =
            (params.recipe_scalar_field_signs[s] < 0) ? -1 : 1;
    }
    return pump;
}

inline void eval_rhs(amrex::MultiFab &a_soln, amrex::MultiFab &a_rhs,
                     const SimulationParameters &params, double dx,
                     const std::array<double, AMREX_SPACEDIM> &center,
                     double time)
{
    const auto &soln_c_arrs = a_soln.const_arrays();
    const auto &rhs_arrs    = a_rhs.arrays();

    if (uses_independent_scalars(params))
    {
        const RLMatterPumpParams pump =
            build_rl_pump(params, params.recipe_num_scalar_fields);

        GRTresnaIndependentScalars matter(
            params.recipe_num_scalar_fields, params.recipe_scalar_field_signs,
            params.recipe_scalar_mass, params.recipe_scalar_lambda, pump);
        CCZ4RHSWithMatter<GRTresnaIndependentScalars,
                          MovingPunctureGaugeWithMatter, FourthOrderDerivatives>
            ccz4rhs(matter, params.ccz4_params, dx, params.sigma,
                    params.formulation, 1.0, center, time);
        amrex::ParallelFor(
            a_rhs, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
            { ccz4rhs(i, j, k, rhs_arrs[box_no], soln_c_arrs[box_no]); });
        return;
    }

    if (uses_complex_scalar(params))
    {
        const RLMatterPumpParams pump = build_rl_pump(params, 1);

        ComplexScalarField matter(params.recipe_scalar_mass,
                                  params.recipe_scalar_lambda,
                                  params.recipe_scalar_sign, pump);
        CCZ4RHSWithMatter<ComplexScalarField,
                          MovingPunctureGaugeWithMatter, FourthOrderDerivatives>
            ccz4rhs(matter, params.ccz4_params, dx, params.sigma,
                    params.formulation, 1.0, center, time);
        amrex::ParallelFor(
            a_rhs, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
            { ccz4rhs(i, j, k, rhs_arrs[box_no], soln_c_arrs[box_no]); });
        return;
    }

    if (uses_bicomplex_scalar(params))
    {
        const RLMatterPumpParams pump = build_rl_pump(params, 2);

        BiComplexScalarField matter(params.recipe_scalar_mass,
                                    params.recipe_scalar_lambda, pump);
        CCZ4RHSWithMatter<BiComplexScalarField,
                          MovingPunctureGaugeWithMatter, FourthOrderDerivatives>
            ccz4rhs(matter, params.ccz4_params, dx, params.sigma,
                    params.formulation, 1.0, center, time);
        amrex::ParallelFor(
            a_rhs, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
            { ccz4rhs(i, j, k, rhs_arrs[box_no], soln_c_arrs[box_no]); });
        return;
    }

    if (params.recipe_exotic_matter)
    {
        ExoticScalarField<DefaultPotential> exotic_scalar(
            DefaultPotential(), params.recipe_support_strength);
        CCZ4RHSWithMatter<ExoticScalarField<DefaultPotential>,
                          MovingPunctureGaugeWithMatter, FourthOrderDerivatives>
            ccz4rhs(exotic_scalar, params.ccz4_params, dx, params.sigma,
                    params.formulation, 1.0, center, time);
        amrex::ParallelFor(
            a_rhs, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
            { ccz4rhs(i, j, k, rhs_arrs[box_no], soln_c_arrs[box_no]); });
        return;
    }

    ScalarField<DefaultPotential> scalar_field;
    CCZ4RHSWithMatter<ScalarField<DefaultPotential>,
                      MovingPunctureGaugeWithMatter, FourthOrderDerivatives>
        ccz4rhs(scalar_field, params.ccz4_params, dx, params.sigma,
                params.formulation, 1.0, center, time);
    amrex::ParallelFor(
        a_rhs, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
        { ccz4rhs(i, j, k, rhs_arrs[box_no], soln_c_arrs[box_no]); });
}

} // namespace RadialRecipeMatter

#endif /* RADIALRECIPE_MATTER_DISPATCH_HPP_ */
