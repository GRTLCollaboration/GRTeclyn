#ifndef RADIALRECIPE_MATTER_DISPATCH_HPP_
#define RADIALRECIPE_MATTER_DISPATCH_HPP_

#include "CCZ4RHSWithMatter.hpp"
#include "ConstraintsWithMatter.hpp"
#include "ExoticScalarField.hpp"
#include "GRTresnaIndependentScalars.hpp"
#include "MovingPunctureGaugeWithMatter.hpp"
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

inline void setup_derived_quantities(int state_index,
                                     const SimulationParameters &params)
{
    if (uses_independent_scalars(params))
    {
        ConstraintsWithMatter<GRTresnaIndependentScalars>::set_up(state_index);
        Weyl4WithMatter<GRTresnaIndependentScalars>::set_up(state_index);
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

inline void eval_rhs(amrex::MultiFab &a_soln, amrex::MultiFab &a_rhs,
                     const SimulationParameters &params, double dx,
                     const std::array<double, AMREX_SPACEDIM> &center,
                     double time)
{
    const auto &soln_c_arrs = a_soln.const_arrays();
    const auto &rhs_arrs    = a_rhs.arrays();

    if (uses_independent_scalars(params))
    {
        GRTresnaIndependentScalars matter(
            params.recipe_num_scalar_fields, params.recipe_scalar_field_signs,
            params.recipe_scalar_mass);
        CCZ4RHSWithMatter<GRTresnaIndependentScalars,
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
