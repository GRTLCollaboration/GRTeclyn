// Catch2 header
#include "catch_amalgamated.hpp"

// Common test headers
#include "InitialData.hpp"

// GRTeclyn headers
#include "CCZ4RHS.hpp"
#include "FourthOrderDerivatives.hpp"

// Old GRTeclyn headers for comparison
#include "CCZ4RHS-fdf5a7a.hpp"
#include "FourthOrderDerivatives-fdf5a7a.hpp"

// AMReX headers
#include "AMReX.H"
#include "AMReX_MultiFab.H"
#include "AMReX_PlotFileUtil.H"

TEST_CASE("CCZ4 RHS")
{
    amrex::Initialize(MPI_COMM_WORLD);
    {
        constexpr int num_cells  = 32;
        constexpr int num_ghosts = 3;
        constexpr amrex::Real dx = 0.5 / (num_cells - 1);

        amrex::Box box(
            amrex::IntVect(0, 0, 0),
            amrex::IntVect(num_cells - 1, num_cells - 1, num_cells - 1));

        amrex::Box ghosted_box = box;
        ghosted_box.grow(num_ghosts);

        amrex::BoxArray box_array{box};
        amrex::DistributionMapping distribution_mapping{box_array};
        amrex::MFInfo mf_info;
        mf_info.SetArena(amrex::The_Managed_Arena());

        amrex::MultiFab in_mf{box_array, distribution_mapping, NUM_CCZ4_VARS,
                              num_ghosts, mf_info};

        // amrex::FArrayBox in_fab{ghosted_box, NUM_CCZ4_VARS,
        //                         amrex::The_Managed_Arena()};

        const auto &in_arrays = in_mf.arrays();
        amrex::ParallelFor(
            in_mf, in_mf.nGrowVect(),
            [=] AMREX_GPU_DEVICE(int ibox, int i, int j, int k)
            {
                const amrex::IntVect iv{i, j, k};
                const amrex::RealVect coords = amrex::RealVect{iv} * dx;

                random_ccz4_initial_data(iv, in_arrays[ibox], coords);
            });

        CCZ4_params_t<MovingPunctureGauge::params_t> current_ccz4_params;
        current_ccz4_params.kappa1            = 0.1;
        current_ccz4_params.kappa2            = 0;
        current_ccz4_params.kappa3            = 1;
        current_ccz4_params.covariantZ4       = true;
        current_ccz4_params.lapse_advec_coeff = 0.0;
        current_ccz4_params.lapse_power       = 1.0;
        current_ccz4_params.lapse_coeff       = 2.0;
        current_ccz4_params.shift_Gamma_coeff = 0.75;
        current_ccz4_params.shift_advec_coeff = 0;
        current_ccz4_params.eta               = 1.82;

        Old::CCZ4_params_t<Old::MovingPunctureGauge::params_t> old_ccz4_params;
        old_ccz4_params.kappa1      = current_ccz4_params.kappa1;
        old_ccz4_params.kappa2      = current_ccz4_params.kappa2;
        old_ccz4_params.kappa3      = current_ccz4_params.kappa3;
        old_ccz4_params.covariantZ4 = current_ccz4_params.covariantZ4;
        old_ccz4_params.lapse_advec_coeff =
            current_ccz4_params.lapse_advec_coeff;
        old_ccz4_params.lapse_power = current_ccz4_params.lapse_power;
        old_ccz4_params.lapse_coeff = current_ccz4_params.lapse_coeff;
        old_ccz4_params.shift_Gamma_coeff =
            current_ccz4_params.shift_Gamma_coeff;
        old_ccz4_params.shift_advec_coeff =
            current_ccz4_params.shift_advec_coeff;
        old_ccz4_params.eta = current_ccz4_params.eta;

        amrex::Real sigma = 0.3;

        CCZ4RHS<MovingPunctureGauge, FourthOrderDerivatives> current_ccz4_rhs{
            current_ccz4_params, dx, sigma};

        Old::CCZ4RHS<Old::MovingPunctureGauge, Old::FourthOrderDerivatives>
            old_ccz4_rhs{old_ccz4_params, dx, sigma};

        // amrex::FArrayBox current_out_fab{box, NUM_CCZ4_VARS,
        //                                  amrex::The_Managed_Arena()};
        // amrex::FArrayBox old_out_fab{box, NUM_CCZ4_VARS,
        //                              amrex::The_Managed_Arena()};
        // amrex::FArrayBox diff_fab{box, NUM_CCZ4_VARS,
        //                           amrex::The_Managed_Arena()};

        amrex::MultiFab current_out_mf{box_array, distribution_mapping,
                                       NUM_CCZ4_VARS, 0, mf_info};
        amrex::MultiFab old_out_mf{box_array, distribution_mapping,
                                   NUM_CCZ4_VARS, 0, mf_info};
        amrex::MultiFab diff_mf{box_array, distribution_mapping, NUM_CCZ4_VARS,
                                0, mf_info};

        const auto &in_c_arrays        = in_mf.const_arrays();
        const auto &current_out_arrays = current_out_mf.arrays();
        const auto &old_out_arrays     = old_out_mf.arrays();
        const auto &diff_arrays        = diff_mf.arrays();

        // Do the current and old CCZ4RHS calculation in the same loop
        amrex::ParallelFor(
            current_out_mf,
            [=] AMREX_GPU_DEVICE(int ibox, int i, int j, int k)
            {
                current_ccz4_rhs.compute(i, j, k, current_out_arrays[ibox],
                                         in_c_arrays[ibox]);
                old_ccz4_rhs.compute(i, j, k, old_out_arrays[ibox],
                                     in_c_arrays[ibox]);

                for (int ivar = 0; ivar < NUM_CCZ4_VARS; ++ivar)
                {
                    diff_arrays[ibox](i, j, k, ivar) =
                        std::fabs(current_out_arrays[ibox](i, j, k, ivar) -
                                  old_out_arrays[ibox](i, j, k, ivar));
                }
            });

        // GPU barrier
        amrex::Gpu::streamSynchronize();

        // write to plot file
        amrex::RealVect dx_Vect{dx, dx, dx};
        amrex::RealBox real_box{box, dx_Vect.dataPtr(),
                                amrex::RealVect::Zero.dataPtr()};

        int coord_sys = 0; // Cartesian
        amrex::Geometry geom{box, &real_box, coord_sys};
        amrex::Vector<std::string> ccz4_var_names_Vect{
            UserVariables::ccz4_variable_names.begin(),
            UserVariables::ccz4_variable_names.end()};

        amrex::WriteSingleLevelPlotfile("pltcur", current_out_mf,
                                        ccz4_var_names_Vect, geom, 0.0, 0);

        amrex::Real max_diff = 0.0;
        amrex::IntVect max_diff_iv{};

        const int cout_precision = Catch::StringMaker<amrex::Real>::precision;
        for (int ivar = 0; ivar < NUM_CCZ4_VARS; ++ivar)
        {
            // diff_fab.maxIndex<amrex::RunOn::Device>(box, max_diff,
            //                                         max_diff_index, ivar);
            max_diff    = diff_mf.max(ivar, 0, true);
            max_diff_iv = diff_mf.maxIndex(ivar);
            int max_diff_boxindex =
                box_array
                    .intersections(amrex::Box{max_diff_iv, max_diff_iv}, true,
                                   0)
                    .front()
                    .first;

            INFO("Max diff for var "
                 << UserVariables::variable_names[ivar] << ": "
                 << std::setprecision(cout_precision) << max_diff << " at "
                 << max_diff_iv);
            INFO("Old value: "
                 << std::setprecision(cout_precision)
                 << old_out_arrays[max_diff_boxindex](max_diff_iv, ivar)
                 << ", Current value: "
                 << current_out_arrays[max_diff_boxindex](max_diff_iv, ivar));
            CHECK_THAT(max_diff, Catch::Matchers::WithinAbs(0.0, 1e-14));
        }

        // GPU barrier
        amrex::Gpu::streamSynchronize();
    }
    amrex::Finalize();
}