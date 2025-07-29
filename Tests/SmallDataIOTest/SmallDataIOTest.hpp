/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */
#ifndef SMALLDATAIOTEST_HPP_
#define SMALLDATAIOTEST_HPP_

// Common includes
#include "doctestCLIArgs.hpp"

// GRTeclyn includes
#include "SmallDataIO.hpp"
#include "SmallDataIOReader.hpp"

// AMReX includes
#include <AMReX.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Gpu.H>
#include <AMReX_Print.H>
#include <AMReX_Random.H>
#include <AMReX_String.H>

// CPP includes
#include <string>
#include <vector>

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
bool check_almost_equal(std::vector<double> vector_1,
                        std::vector<double> vector_2, double err_tol);

std::vector<double> generate_random_numbers(const int Npts);

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
void test_small_data_io_writer(
    const std::vector<SmallDataIOReader::column_t> &col,
    const int data_precision);

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
std::vector<SmallDataIOReader::column_t>
test_small_data_io_reader(const std::vector<std::string> &column_names);

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
std::vector<SmallDataIOReader::column_t>
test_small_data_io_reader(const int a_min_col, const int a_max_col);

std::vector<double> test_small_data_io_reader(const int a_col);

void run_small_data_io_test();

#endif /* SMALLDATAIOTEST_HPP_ */
