#ifndef RADIALRECIPE_MATTER_DISPATCH_HPP_
#define RADIALRECIPE_MATTER_DISPATCH_HPP_

#include "BiComplexScalarField.hpp"
#include "CCZ4RHSWithMatter.hpp"
#include "ComplexScalarField.hpp"
#include "ConstraintsWithMatter.hpp"
#include "ControllerReservoirMatter.hpp"
#include "ExoticScalarField.hpp"
#include "GRTresnaIndependentScalars.hpp"
#include "MovingPunctureGaugeWithMatter.hpp"
#include "RadialRecipeMatterConstraints.hpp"
#include "RLMatterPumpParams.hpp"
#include "ScalarField.hpp"
#include "SimulationParameters.hpp"
#include "SpongeZone.hpp"
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
                                        int num_fields, double time = 0.0)
{
    RLMatterPumpParams pump;
    int n = params.rl_num_lumps;
    if (n > GRTRESNA_MAX_INDEPENDENT_SCALARS)
        n = GRTRESNA_MAX_INDEPENDENT_SCALARS;
    if (n < 0)
        n = 0;
    // Transient igniter: after rl_pump_stop_time the pump is fully off so the
    // measured window is a conservative Einstein-Klein-Gordon evolution.
    const bool pump_stopped =
        (params.rl_pump_stop_time >= 0.0 && time >= params.rl_pump_stop_time);
    if (pump_stopped)
        n = 0;
    pump.num_sites       = n;
    pump.width           = params.rl_pump_width;
    pump.governor_center = params.rl_l2_ham_governor_center;
    pump.governor_width  = params.rl_l2_ham_governor_width;
    pump.governor        = RLRuntime::tanh_governor(
        RLRuntime::cached_L2_Ham(), pump.governor_center, pump.governor_width);
    pump.num_fields      = num_fields;
    // Closed-loop PD "trap" controller gains (0 => legacy open-loop source).
    pump.k_p             = pump_stopped ? 0.0 : params.rl_pump_kp;
    pump.k_d             = pump_stopped ? 0.0 : params.rl_pump_kd;
    pump.target_profile  = params.rl_pump_target_profile;
    pump.target_width    = params.rl_pump_target_width;
    pump.target_amp      = params.rl_pump_target_amp;
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

//! Fill constraints for the active matter model. When
//! controller_reservoir_mode >= 1 the reservoir is included in the diagnostic
//! EMT (ledger / backreaction modes). ``with_abs_terms`` runs a second pass
//! writing Ham_abs / Mom_abs into comps 4 / 5-7 (cst must have >= 8 comps).
template <class Matter>
inline void fill_constraints_maybe_abs(
    amrex::MultiFab &cst, const amrex::MultiFab &state_new, const Matter &matter,
    amrex::Real dx0, const std::array<double, AMREX_SPACEDIM> &center,
    amrex::Real time, bool with_abs_terms)
{
    fill_matter_constraints(cst, state_new, matter, dx0, center, time);
    if (with_abs_terms)
    {
        fill_matter_constraint_abs_terms(cst, state_new, matter, dx0, center,
                                         time);
    }
}

//! ``reservoir_mode_override >= 0`` forces a specific reservoir mode instead of
//! reading it from the params.  This exists so the SAFETY GOVERNOR can be fed
//! the true, physical constraint norm.  Feeding it the ledger-corrected value
//! couples the controller's own bookkeeping back into its safety interlock: if
//! the ledger misbehaves, the governor throttles the pump on the strength of a
//! number that describes the ledger rather than the spacetime.  That is exactly
//! what happened in the first always-on campaign (the reservoir diverged, the
//! reported Ham crossed the 0.035 governor threshold at t~7.5, and the pump was
//! cut to zero by t~10 in every mode-1 run).  Governor input must stay physical.
inline void fill_active_constraints(
    amrex::MultiFab &cst, const amrex::MultiFab &state_new,
    const SimulationParameters &params, amrex::Real dx0,
    const std::array<double, AMREX_SPACEDIM> &center, amrex::Real time,
    bool with_abs_terms = false, int reservoir_mode_override = -1)
{
    const int mode = (reservoir_mode_override >= 0)
                         ? reservoir_mode_override
                         : params.controller_reservoir_mode;
    // Constraints path: include reservoir in EMT for both ledger (1) and
    // backreaction (2) so the reported norms are pump-corrected.
    const bool include_em = (mode >= 1);

    if (uses_independent_scalars(params))
    {
        GRTresnaIndependentScalars matter(
            params.recipe_num_scalar_fields, params.recipe_scalar_field_signs,
            params.recipe_scalar_mass, params.recipe_scalar_lambda);
        fill_constraints_maybe_abs(cst, state_new, matter, dx0, center, time,
                                   with_abs_terms);
    }
    else if (uses_complex_scalar(params))
    {
        const RLMatterPumpParams pump = build_rl_pump(params, 1, time);
        ComplexScalarField inner(params.recipe_scalar_mass,
                                 params.recipe_scalar_lambda,
                                 params.recipe_scalar_sign,
                                 params.recipe_scalar_mu, pump);
        if (mode >= 1)
        {
            using Wrapped =
                ControllerReservoirMatter<ComplexScalarField, false>;
            Wrapped matter(inner, pump, mode, include_em);
            fill_constraints_maybe_abs(cst, state_new, matter, dx0, center,
                                       time, with_abs_terms);
        }
        else
        {
            fill_constraints_maybe_abs(cst, state_new, inner, dx0, center, time,
                                       with_abs_terms);
        }
    }
    else if (uses_bicomplex_scalar(params))
    {
        const RLMatterPumpParams pump = build_rl_pump(params, 2, time);
        BiComplexScalarField inner(params.recipe_scalar_mass,
                                   params.recipe_scalar_lambda,
                                   params.recipe_scalar_mu, pump);
        if (mode >= 1)
        {
            using Wrapped =
                ControllerReservoirMatter<BiComplexScalarField, true>;
            Wrapped matter(inner, pump, mode, include_em);
            fill_constraints_maybe_abs(cst, state_new, matter, dx0, center,
                                       time, with_abs_terms);
        }
        else
        {
            fill_constraints_maybe_abs(cst, state_new, inner, dx0, center, time,
                                       with_abs_terms);
        }
    }
    else if (params.recipe_exotic_matter)
    {
        ExoticScalarField<DefaultPotential> exotic_scalar(
            DefaultPotential(), params.recipe_support_strength);
        fill_constraints_maybe_abs(cst, state_new, exotic_scalar, dx0, center,
                                   time, with_abs_terms);
    }
    else
    {
        ScalarField<DefaultPotential> scalar_field;
        fill_constraints_maybe_abs(cst, state_new, scalar_field, dx0, center,
                                   time, with_abs_terms);
    }
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
            build_rl_pump(params, params.recipe_num_scalar_fields, time);

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
    }
    else if (uses_complex_scalar(params))
    {
        const RLMatterPumpParams pump = build_rl_pump(params, 1, time);
        const int mode = params.controller_reservoir_mode;
        // Mode 2: include reservoir in the Einstein RHS (backreaction).
        // Mode 1: evolve reservoir but do not source the metric from it.
        const bool include_em = (mode >= 2);

        ComplexScalarField inner(params.recipe_scalar_mass,
                                 params.recipe_scalar_lambda,
                                 params.recipe_scalar_sign,
                                 params.recipe_scalar_mu, pump);
        if (mode >= 1)
        {
            using Wrapped =
                ControllerReservoirMatter<ComplexScalarField, false>;
            Wrapped matter(inner, pump, mode, include_em);
            CCZ4RHSWithMatter<Wrapped, MovingPunctureGaugeWithMatter,
                              FourthOrderDerivatives>
                ccz4rhs(matter, params.ccz4_params, dx, params.sigma,
                        params.formulation, 1.0, center, time);
            amrex::ParallelFor(
                a_rhs, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                { ccz4rhs(i, j, k, rhs_arrs[box_no], soln_c_arrs[box_no]); });
        }
        else
        {
            CCZ4RHSWithMatter<ComplexScalarField, MovingPunctureGaugeWithMatter,
                              FourthOrderDerivatives>
                ccz4rhs(inner, params.ccz4_params, dx, params.sigma,
                        params.formulation, 1.0, center, time);
            amrex::ParallelFor(
                a_rhs, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                { ccz4rhs(i, j, k, rhs_arrs[box_no], soln_c_arrs[box_no]); });
        }
    }
    else if (uses_bicomplex_scalar(params))
    {
        const RLMatterPumpParams pump = build_rl_pump(params, 2, time);
        const int mode = params.controller_reservoir_mode;
        const bool include_em = (mode >= 2);

        BiComplexScalarField inner(params.recipe_scalar_mass,
                                   params.recipe_scalar_lambda,
                                   params.recipe_scalar_mu, pump);
        if (mode >= 1)
        {
            using Wrapped =
                ControllerReservoirMatter<BiComplexScalarField, true>;
            Wrapped matter(inner, pump, mode, include_em);
            CCZ4RHSWithMatter<Wrapped, MovingPunctureGaugeWithMatter,
                              FourthOrderDerivatives>
                ccz4rhs(matter, params.ccz4_params, dx, params.sigma,
                        params.formulation, 1.0, center, time);
            amrex::ParallelFor(
                a_rhs, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                { ccz4rhs(i, j, k, rhs_arrs[box_no], soln_c_arrs[box_no]); });
        }
        else
        {
            CCZ4RHSWithMatter<BiComplexScalarField,
                              MovingPunctureGaugeWithMatter,
                              FourthOrderDerivatives>
                ccz4rhs(inner, params.ccz4_params, dx, params.sigma,
                        params.formulation, 1.0, center, time);
            amrex::ParallelFor(
                a_rhs, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                { ccz4rhs(i, j, k, rhs_arrs[box_no], soln_c_arrs[box_no]); });
        }
    }
    else if (params.recipe_exotic_matter)
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
    }
    else
    {
        ScalarField<DefaultPotential> scalar_field;
        CCZ4RHSWithMatter<ScalarField<DefaultPotential>,
                          MovingPunctureGaugeWithMatter, FourthOrderDerivatives>
            ccz4rhs(scalar_field, params.ccz4_params, dx, params.sigma,
                    params.formulation, 1.0, center, time);
        amrex::ParallelFor(
            a_rhs, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
            { ccz4rhs(i, j, k, rhs_arrs[box_no], soln_c_arrs[box_no]); });
    }

    // Sponge zone: apply extra radially-ramped KO dissipation in the outer
    // shell to absorb outgoing waves before they hit the Sommerfeld boundary.
    if (params.sponge_params.enabled)
    {
        const SpongeZone sponge(params.sponge_params, dx);
        amrex::ParallelFor(
            a_rhs, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
            { sponge.apply(i, j, k, rhs_arrs[box_no], soln_c_arrs[box_no]); });
    }
}

} // namespace RadialRecipeMatter

#endif /* RADIALRECIPE_MATTER_DISPATCH_HPP_ */
