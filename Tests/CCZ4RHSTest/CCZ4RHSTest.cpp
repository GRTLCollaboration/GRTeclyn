// Catch2 header
#include "catch_amalgamated.hpp"

// Common test headers
#include "InitialData.hpp"

// GRAMReX headers
#include "CCZ4RHS.hpp"
#include "FourthOrderDerivatives.hpp"

// Old GRAMReX headers for comparison
#include "CCZ4RHS-fdf5a7a.hpp"
#include "FourthOrderDerivatives-fdf5a7a.hpp"

// AMReX headers
#include "AMReX.H"
#include "AMReX_FArrayBox.H"

TEST_CASE("CCZ4 RHS")
{
    amrex::Initialize(MPI_COMM_WORLD);
    constexpr int num_cells  = 32;
    constexpr int num_ghosts = 3;
    constexpr double dx      = 0.5 / (num_cells - 1);

    amrex::Box box(amrex::IntVect(0, 0, 0),
                   amrex::IntVect(num_cells - 1, num_cells - 1, num_cells - 1));

    amrex::Box ghosted_box = box;
    ghosted_box.grow(num_ghosts);

    amrex::FArrayBox in_fab{ghosted_box, NUM_CCZ4_VARS,
                            amrex::The_Managed_Arena()};

    const amrex::Array4<amrex::Real> &in_array = in_fab.array();
    amrex::ParallelFor(ghosted_box,
                       [=] AMREX_GPU_DEVICE(int i, int j, int k)
                       {
                           const amrex::IntVect iv{i, j, k};
                           const amrex::RealVect coords =
                               amrex::RealVect{iv} * dx;

                           random_ccz4_initial_data(iv, in_array, coords);
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
    old_ccz4_params.kappa1            = current_ccz4_params.kappa1;
    old_ccz4_params.kappa2            = current_ccz4_params.kappa2;
    old_ccz4_params.kappa3            = current_ccz4_params.kappa3;
    old_ccz4_params.covariantZ4       = current_ccz4_params.covariantZ4;
    old_ccz4_params.lapse_advec_coeff = current_ccz4_params.lapse_advec_coeff;
    old_ccz4_params.lapse_power       = current_ccz4_params.lapse_power;
    old_ccz4_params.lapse_coeff       = current_ccz4_params.lapse_coeff;
    old_ccz4_params.shift_Gamma_coeff = current_ccz4_params.shift_Gamma_coeff;
    old_ccz4_params.shift_advec_coeff = current_ccz4_params.shift_advec_coeff;
    old_ccz4_params.eta               = current_ccz4_params.eta;

    amrex::Real sigma = 0.3;

    CCZ4RHS<MovingPunctureGauge, FourthOrderDerivatives> current_ccz4_rhs{
        current_ccz4_params, dx, sigma};

    Old::CCZ4RHS<Old::MovingPunctureGauge, Old::FourthOrderDerivatives>
        old_ccz4_rhs{old_ccz4_params, dx, sigma};

    amrex::FArrayBox current_out_fab{box, NUM_CCZ4_VARS,
                                     amrex::The_Managed_Arena()};
    amrex::FArrayBox old_out_fab{box, NUM_CCZ4_VARS,
                                 amrex::The_Managed_Arena()};
    amrex::FArrayBox diff_fab{box, NUM_CCZ4_VARS, amrex::The_Managed_Arena()};

    const auto &in_c_array        = in_fab.const_array();
    const auto &current_out_array = current_out_fab.array();
    const auto &old_out_array     = old_out_fab.array();
    const auto &diff_array        = diff_fab.array();

    // Do the current and old CCZ4RHS calculation in the same loop
    amrex::ParallelFor(
        box,
        [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            current_ccz4_rhs.compute(i, j, k, current_out_array, in_c_array);
            old_ccz4_rhs.compute(i, j, k, old_out_array, in_c_array);

            for (int ivar = 0; ivar < NUM_CCZ4_VARS; ++ivar)
            {
                diff_array(i, j, k, ivar) =
                    std::fabs(current_out_array(i, j, k, ivar) -
                              old_out_array(i, j, k, ivar));
            }
        });

    // GPU barrier
    amrex::Gpu::streamSynchronize();

    amrex::Real max_diff = 0.0;
    amrex::IntVect max_diff_index{};

    const int cout_precision = Catch::StringMaker<amrex::Real>::precision;
    for (int ivar = 0; ivar < NUM_CCZ4_VARS; ++ivar)
    {
        diff_fab.maxIndex(box, max_diff, max_diff_index, ivar);

        INFO("Max diff for var " << UserVariables::variable_names[ivar] << ": "
                                 << std::setprecision(cout_precision)
                                 << max_diff << " at " << max_diff_index);
        INFO("Old value: " << std::setprecision(cout_precision)
                           << old_out_array(max_diff_index, ivar)
                           << ", Current value: "
                           << current_out_array(max_diff_index, ivar));
        CHECK_THAT(max_diff, Catch::Matchers::WithinAbs(0.0, 1e-14));
    }

    amrex::Finalize();
}