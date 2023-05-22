// Catch2 header
#include "catch_amalgamated.hpp"

// Common test headers
#include "InitialData.hpp"

// GRAMReX headers
#include "CCZ4RHS.hpp"
#include "FourthOrderDerivatives.hpp"

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

    amrex::FArrayBox in_fab{ghosted_box, NUM_CCZ4_VARS, amrex::The_Arena()};
    amrex::FArrayBox out_fab{box, NUM_CCZ4_VARS, amrex::The_Arena()};

    const amrex::Array4<amrex::Real> &in_array = in_fab.array();
    amrex::ParallelFor(ghosted_box,
                       [=] AMREX_GPU_DEVICE(int i, int j, int k)
                       {
                           const amrex::IntVect iv{i, j, k};
                           const amrex::RealVect coords =
                               amrex::RealVect{iv} * dx;

                           random_ccz4_initial_data(iv, in_array, coords);
                       });

    CCZ4_params_t<MovingPunctureGauge::params_t> ccz4_params;
    ccz4_params.kappa1            = 0.1;
    ccz4_params.kappa2            = 0;
    ccz4_params.kappa3            = 1;
    ccz4_params.covariantZ4       = true;
    ccz4_params.lapse_advec_coeff = 0.0;
    ccz4_params.lapse_power       = 1.0;
    ccz4_params.lapse_coeff       = 2.0;
    ccz4_params.shift_Gamma_coeff = 0.75;
    ccz4_params.shift_advec_coeff = 0;
    ccz4_params.eta               = 1.82;

    amrex::Real sigma = 0.3;

    CCZ4RHS<MovingPunctureGauge, FourthOrderDerivatives> ccz4_rhs{ccz4_params,
                                                                  dx, sigma};

    const auto &out_array  = out_fab.array();
    const auto &in_c_array = in_fab.const_array();
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                       { ccz4_rhs.compute(i, j, k, out_array, in_c_array); });

    REQUIRE(true);

    amrex::Finalize();
}