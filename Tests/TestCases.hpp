/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */
#ifndef TESTCASES_HPP_
#define TESTCASES_HPP_

// doctest header
#include "doctest.h"

// AMReX includes
#include <AMReX.H>

// Test cases
#include "BSSNMatterTest.hpp"
#include "CCZ4GeometryUnitTest.hpp"
#include "CCZ4RHSTest.hpp"
#include "ConstraintsTest.hpp"
#include "CoordinateTransformationsTest.hpp"
#include "DerivativeUnitTests.hpp"
#include "EMTensorTest.hpp"
#include "ParticleInterpolatorUnitTest.hpp"
#include "PositiveChiAndLapseUnitTest.hpp"
#include "PunctureTrackerTest.hpp"
#include "SmallDataIOTest.hpp"
#include "SphericalHarmonicTest.hpp"
#include "Weyl4Test.hpp"
#include "Weyl4WithMatterTest.hpp"

TEST_CASE("CCZ4Geometry") { run_ccz4_geometry_unit_tests(); }

TEST_CASE("BSSNMatter"
#ifndef AMREX_USE_HDF5
          * doctest::skip()
#endif
)
{
    run_bssn_matter_test();
}

TEST_CASE("CCZ4 Geometry") { run_ccz4_geometry_unit_tests(); }

TEST_CASE("Particle Interpolator") { run_particle_interpolator_test(); }

TEST_CASE("CCZ4RHS") { run_ccz4_rhs_test(); }

TEST_CASE("Constraints"
#ifndef AMREX_USE_HDF5
          * doctest::skip()
#endif
)
{
    run_constraints_test();
}

TEST_CASE("CoordinateTransformations")
{
    run_coordinate_transformations_test();
}

TEST_CASE("DerivativeUnitTests") { run_derivative_unit_tests(); }

TEST_CASE("EMTensor"
#ifndef AMREX_USE_HDF5
          * doctest::skip()
#endif
)
{
    run_emtensor_test();
}

TEST_CASE("PositiveChiAndLapse") { run_positive_chi_and_lapse_unit_test(); }

TEST_CASE("PunctureTracker") { run_puncture_tracker_test(); }

TEST_CASE("SmallDataIO") { run_small_data_io_test(); }

TEST_CASE("SphericalHarmonics") { run_spherical_harmonic_test(); }

TEST_CASE("Weyl4"
#ifndef AMREX_USE_HDF5
          * doctest::skip()
#endif
)
{
    run_weyl4_test();
}

TEST_CASE("Weyl4WithMatter"
#ifndef AMREX_USE_HDF5
          * doctest::skip()
#endif
)
{
    run_matter_weyl4_test();
}

#endif /* TESTCASES_HPP_ */
