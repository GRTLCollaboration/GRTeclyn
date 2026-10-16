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

#include "AMReX_GpuParallelReduce.H"
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

                state_arrays[ibox](iv, c_phi_Re)    = 0.0;
                state_arrays[ibox](iv, c_phi_Im)    = 0.0;
                state_arrays[ibox](iv, c_polystate) = 0.0;
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
                           amrex::Real a_time) override
    {
        a_rhs.setVal(1.);
    };

    void specific_update_ode(amrex::MultiFab &a_soln) override {};

    void specific_post_timestep() override
    { // Set up the Spectral Wishes

        amrex::MultiFab &state   = get_new_data(state_index);
        amrex::Real current_time = get_state_data(state_index).curTime();
        auto *gr_amr_ptr         = get_gr_amr_ptr();
        SpectralWishes<DefaultScalarField> my_spectral_wishes(
            state_index, current_time, level, state.nGrow());
        const auto spectral_wishes_results =
            my_spectral_wishes.compute(gr_amr_ptr);
        const auto &wishes = my_spectral_wishes.get_wishes();
        const auto labels  = my_spectral_wishes.get_labels();

        // Make sure the parameter file actually asks for something to test
        CHECK(!wishes.empty());
        CHECK(spectral_wishes_results.size() == wishes.size());

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

        const auto num_cells = static_cast<amrex::Real>(countCells());

        // Mean of every component of derive_mf
        amrex::Vector<amrex::Real> mean_direct_calc = sum_powers(
            derive_mf, amrex::Vector<amrex::Real>(num_comp_total), 1);
        for (auto &mean : mean_direct_calc)
        {
            mean /= num_cells;
        }

        // Variance of every component of derive_mf
        amrex::Vector<amrex::Real> variance_direct_calc =
            sum_powers(derive_mf, mean_direct_calc, 2);
        for (auto &variance : variance_direct_calc)
        {
            variance /= num_cells - 1.0;
        }

        // Compare every wish with the direct calculation of the same
        // quantity
        for (amrex::Long iwish = 0; iwish < wishes.size(); ++iwish)
        {
            const auto &wish = wishes[iwish];
            const int comp   = derive_mf_comp(wish.variable, dcomp);

            amrex::Real direct_calc{0.};
            switch (wish.operation)
            {
            case SpectralWishOperation::Mean:
                direct_calc = mean_direct_calc[comp];
                break;
            case SpectralWishOperation::Variance:
                direct_calc = variance_direct_calc[comp];
                break;
            }

            amrex::Print() << labels[iwish] << " (spectral wishes): "
                           << spectral_wishes_results[iwish] << "\n";
            amrex::Print() << labels[iwish]
                           << " (direct calculation): " << direct_calc << "\n";

            INFO("Spectral wish: ", labels[iwish]);
            CHECK(spectral_wishes_results[iwish] ==
                  doctest::Approx(direct_calc).epsilon(1e-12));
        }
    };

    //! Position of a variable in derive_mf: the state variables come first,
    //! followed by the constraints starting at a_constraints_comp
    static int derive_mf_comp(const std::string &a_variable,
                              int a_constraints_comp)
    {
        int state_type{0};
        int state_comp{0};
        if (amrex::AmrLevel::isStateVariable(a_variable, state_type,
                                             state_comp))
        {
            REQUIRE(state_comp < a_constraints_comp);
            return state_comp;
        }

        const auto *derive_rec =
            amrex::AmrLevel::get_derive_lst().get(a_variable);
        REQUIRE(derive_rec != nullptr);
        for (int icomp = 0; icomp < derive_rec->numDerive(); ++icomp)
        {
            if (derive_rec->variableName(icomp) == a_variable)
            {
                return a_constraints_comp + icomp;
            }
        }

        FAIL("Variable ", a_variable, " not found in derive_mf");
        return -1;
    }

    //! For every component, the sum over all cells (and MPI ranks) of
    //! (value - a_shift[component])^a_power
    static amrex::Vector<amrex::Real>
    sum_powers(const amrex::MultiFab &a_mf,
               const amrex::Vector<amrex::Real> &a_shift, int a_power)
    {
        const int num_comps = a_mf.nComp();

        // Use a GPU safe data container to transfer data from host to device
        // and vice versa
        amrex::Gpu::DeviceVector<amrex::Real> device_sums(num_comps, 0.0);
        amrex::Gpu::DeviceVector<amrex::Real> device_shift(num_comps);
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, a_shift.begin(),
                         a_shift.end(), device_shift.begin());

        auto *device_sums_ptr        = device_sums.dataPtr();
        const auto *device_shift_ptr = device_shift.dataPtr();

        const auto &mf_arrays = a_mf.const_arrays();

        // Sum over all cell values of a_mf
        amrex::ParallelFor(
            a_mf,
            [=] AMREX_GPU_DEVICE(int ibox, int ix, int iy, int iz)
            {
                for (int component = 0; component < num_comps; ++component)
                {
                    const amrex::Real deviation =
                        mf_arrays[ibox](ix, iy, iz, component) -
                        device_shift_ptr[component];
                    amrex::Real term{1.};
                    for (int ipower = 0; ipower < a_power; ++ipower)
                    {
                        term *= deviation;
                    }
                    // Unlike Gpu::Atomic this is also atomic on the host
                    // when using OpenMP
                    amrex::HostDevice::Atomic::Add(&device_sums_ptr[component],
                                                   term);
                }
            });
        // GPU barrier
        amrex::Gpu::streamSynchronize();

        // Copy back from device (ParallelFor kernel) to host
        // This is blocking - ok because the result is immediately being used
        // otherwise use Gpu::AsyncCopy
        amrex::Vector<amrex::Real> sums(num_comps);
        amrex::Gpu::copy(amrex::Gpu::deviceToHost, device_sums.begin(),
                         device_sums.end(), sums.begin());

        // If running with MPI each rank will hold a different portion of
        // a_mf so need to reduce (sum) into a single value across all grids
        amrex::ParallelAllReduce::Sum(
            sums.data(), num_comps, amrex::ParallelDescriptor::Communicator());

        return sums;
    }

    void tag_cells(amrex::TagBoxArray &a_tag_box_array,
                   amrex::Real a_regrid_threshold) final {};
};

#endif /* SPECTRALWISHESTESTLEVEL_HPP_ */
