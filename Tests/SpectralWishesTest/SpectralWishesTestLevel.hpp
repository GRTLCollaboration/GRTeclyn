/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SPECTRALWISHESTESTLEVEL_HPP_
#define SPECTRALWISHESTESTLEVEL_HPP_

#include "CCZ4RHSWithMatter.hpp"
#include "ConstraintsWithMatter.hpp"
#include "DefaultLevelBld.hpp"
#include "DefaultPotential.hpp"
#include "EMTensor.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRAmrLevel.hpp"
#include "GRParmParse.hpp"
#include "MovingPunctureGauge.hpp"
#include "ScalarField.hpp"

#include "InitialData.hpp"
#include "SpectralWishes.hpp"

#include "AMReX_ParallelDescriptor.H"

/// Evolution level for a real scalar field minimally coupled to gravity.
class SpectralWishesTestLevel : public GRAmrLevel
{
  public:
    using GRAmrLevel::GRAmrLevel;

    using DefaultScalarField =
        ScalarField<DefaultPotential, FourthOrderDerivatives>;

    using ScalarFieldConstraints = ConstraintsWithMatter<DefaultScalarField>;

    static void variableSetUp()
    {
        state_variable_set_up();
        ScalarFieldConstraints::set_up(state_index);
    };

    void specific_advance() override {};

    void initData() override
    {
        GRParmParse pp;

        auto geometry            = Geom();
        const auto dx            = geometry.CellSizeArray();
        amrex::MultiFab &state   = get_new_data(state_index);
        const auto &state_arrays = state.arrays();
        // NOLINTBEGIN(bugprone-easily-swappable-parameters)
        amrex::ParallelFor(
            state, state.nGrowVect(),
            [=] AMREX_GPU_DEVICE(int ibox, int i, int j, int k)
            // NOLINTEND(bugprone-easily-swappable-parameters)
            {
                const amrex::IntVect iv{i, j, k};
                const amrex::RealVect coords = amrex::RealVect{iv} * dx[0];
                amrex::Real x                = coords[0];
                amrex::Real y                = coords[1];
                amrex::Real z                = coords[2];

                random_ccz4_initial_data(iv, state_arrays[ibox], coords);

                random_matter_bssn_initial_data(iv, state_arrays[ibox], coords);
            });

        pp.add("ccz4.kappa1", 0.0);
        pp.add("ccz4.kappa2", 0.0);
        pp.add("ccz4.kappa3", 0.0);
        pp.add("ccz4.covariantZ4", true);

        pp.add("gauge.shift_Gamma_coeff", 0.75);
        pp.add("gauge.lapse_advec_coeff", 1.0);
        pp.add("gauge.lapse_power", 1.0);
        pp.add("gauge.lapse_coeff", 2.0);
        pp.add("gauge.shift_advec_coeff", 0.0);
        pp.add("gauge.eta", 1.0);

        pp.add("evolution.sigma", 0.1);
        pp.add("ccz4.formulation", CCZ4RHS<>::USE_BSSN);
    };

    void specific_eval_rhs(amrex::MultiFab &a_soln, amrex::MultiFab &a_rhs,
                           amrex::Real a_time) override {};

    void specific_update_ode(amrex::MultiFab &a_soln) override {};

    void specific_post_timestep() override
    { // Set up the Spectral Wishes

        amrex::MultiFab &state   = get_new_data(state_index);
        const auto &state_arrays = state.arrays();

        amrex::Real current_time = get_state_data(state_index).curTime();
        auto *gr_amr_ptr         = get_gr_amr_ptr();
        SpectralWishes<DefaultScalarField> my_spectral_wishes(
            state_index, current_time, level, state.nGrow());
        const amrex::Real mean_spectral_wishes =
            my_spectral_wishes.compute(gr_amr_ptr);

        // Calculate the Wishes directly from the MultiFabs for comparison

        int *bcrec = nullptr;

        // This is not the same as NUM_VARS ( = 30)
        constexpr int num_bssn_matter_vars = c_Pi + 1;
        constexpr int dcomp                = num_bssn_matter_vars;

        constexpr int num_comp_constraints = 1 + AMREX_SPACEDIM; // ham + moms;
        constexpr int num_comp_total =
            num_bssn_matter_vars + num_comp_constraints;

        amrex::MultiFab derive_mf{state.boxArray(), state.DistributionMap(),
                                  num_comp_total, 0};
        derive_mf.setVal(0.0); // initialise to zero

        // Copy over the values from state multifab so ParallelFor only has
        // to be launched once for both state and derived variables
        amrex::MultiFab::Copy(derive_mf, state, 0, 0, num_bssn_matter_vars, 0);
        ConstraintsWithMatter<DefaultScalarField>::compute_mf(
            derive_mf, dcomp, num_comp_constraints, state, geom, current_time,
            bcrec, level);

        amrex::Vector<amrex::Real> mean_direct_calc(derive_mf.nComp());
        std::ranges::fill(mean_direct_calc.begin(), mean_direct_calc.end(),
                          0.0);

        // Use AsyncArrays to transfer data from host to device and vice versa
        amrex::AsyncArray<amrex::Real> async_arr(mean_direct_calc.data(),
                                                 derive_mf.nComp());
        auto *async_arr_data_ptr =
            async_arr.data(); // data pointer to data contained in async_arr

        const auto &derive_mf_arrays = derive_mf.arrays();
        amrex::ParallelFor(
            derive_mf,
            [=] AMREX_GPU_DEVICE(int ibox, int ix, int iy, int iz)
            {
                const amrex::CellData<amrex::Real> cell =
                    derive_mf_arrays[ibox].cellData(ix, iy, iz);
                for (int component = 0; component < num_comp_total; ++component)
                {
                    async_arr_data_ptr[component] += cell[component];
                }
            });

        // Copy back from device (ParallelFor kernel) to host
        async_arr.copyToHost(mean_direct_calc.data(), derive_mf.nComp());

        for (int component = 0; component < derive_mf.nComp(); ++component)
        {
            mean_direct_calc[component] /=
                static_cast<amrex::Real>(countCells());
            amrex::ParallelDescriptor::ReduceRealSum(
                mean_direct_calc[component]);
        }

        amrex::Print() << "Mean (direct calculation): " << mean_direct_calc[29]
                       << "\n";

        // GPU barrier
        amrex::Gpu::streamSynchronize();

        CHECK(mean_direct_calc[29] ==
              doctest::Approx(mean_spectral_wishes).epsilon(1e-12));
    };

    void tag_cells(amrex::TagBoxArray &a_tag_box_array,
                   amrex::Real a_regrid_threshold) final {};
};

#endif /* SPECTRALWISHESTESTLEVEL_HPP_ */
