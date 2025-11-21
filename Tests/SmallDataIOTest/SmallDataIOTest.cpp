/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// Doctest header
#include "doctest.h"

// Test header
#include "SmallDataIOTest.hpp"

std::vector<double> generate_random_numbers(const int Npts)
{
    amrex::Gpu::DeviceVector<double> data_device(Npts);
    amrex::FillRandom(data_device.begin(), Npts);
    std::vector<double> data_host(Npts);

    amrex::Gpu::copyAsync(amrex::Gpu::deviceToHost, data_device.begin(),
                          data_device.end(), data_host.begin());

    // NB: Remember to call streamSynchronize once all copyAsyncs are done

    return data_host;
}
// NOLINTBEGIN(bugprone-easily-swappable-parameters)
bool check_almost_equal(std::vector<double> vector_1,
                        std::vector<double> vector_2, double err_tol)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    // Check that the vector_1 is equal to vector_2 within err_tol

    // Note that we couldn't have used amrex::almostEqual here because that
    // assumes that the two values are floats not vectors

    AMREX_ASSERT_WITH_MESSAGE(
        vector_1.size() == vector_2.size(),
        "Vectors must be the same length for the comparison to work!\n");

    for (int i = 0; i < vector_1.size(); ++i)
    {
        CHECK_MESSAGE(vector_1[i] ==
                          doctest::Approx(vector_2[i]).epsilon(err_tol),
                      "i= " << i);
    }

    return true;
}
// NOLINTBEGIN(bugprone-easily-swappable-parameters)
void test_small_data_io_writer(const std::vector<SmallDataIO::column_t> &col,
                               const int data_precision)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    const std::string filename_prefix{"test_"};
    amrex::Real dt{0.1};
    amrex::Real time{0.0};
    amrex::Real restart_time{0.0};
    bool first_step{true};

    SmallDataIO test_file(filename_prefix, dt, time, restart_time,
                          SmallDataIO::NEW, first_step, ".dat", data_precision);

    static const std::vector<std::string> header1_strings{"x", "y", "z"};

    // write a header containing labels for each column
    test_file.write_header_line(header1_strings); // pre-header string is time

    for (int i = 0; i < col[0].size(); ++i)
    {
        std::vector<amrex::Real> some_data{col[0][i], col[1][i], col[2][i]};
        test_file.write_time_data_line(some_data);
    }
}

std::vector<SmallDataIO::column_t>
test_small_data_io_reader(const std::vector<std::string> &column_names)
{
    SmallDataIO test_reader("test_000000");

    // Could print out file structure as well
    // test_reader.print_file_structure();

    std::vector<std::string> header;
    test_reader.get_header_strings(header, 0);

    // If a column name doesn't exist in the header, amrex::Abort will
    // be called

    std::vector<SmallDataIO::column_t> data;
    test_reader.get_columns(data, column_names, 0);

    SmallDataIO::broadcast_data(data);

    return data;
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
std::vector<SmallDataIO::column_t>
test_small_data_io_reader(const int a_min_col, const int a_max_col)
{
    // If we know the data structure already we can pass it in like so:

    const SmallDataIO::file_structure_t test_file_structure{
        1, {0}, {1}, {100}, {4}};

    SmallDataIO test_reader("test_000000", &test_file_structure);

    std::vector<SmallDataIO::column_t> data;
    test_reader.get_columns(data, a_min_col, a_max_col, 0);

    SmallDataIO::broadcast_data(data);

    return data;
}

std::vector<double> test_small_data_io_reader(const int a_col)

{
    SmallDataIO test_reader("test_000000");

    std::vector<SmallDataIO::column_t> data;
    test_reader.get_column(data, a_col, 0);

    return data[0];
}

void run_small_data_io_test()
{
    int amrex_argc    = doctest::cli_args.argc();
    char **amrex_argv = doctest::cli_args.argv();

    // NOLINTNEXTLINE(bugprone-casting-through-void) // Open MPI triggers this
    amrex::Initialize(amrex_argc, amrex_argv);
    {

        // Generate some random numbers to write out
        static const int Npts{100};

        // This is GPU safe, but must streamSynchronize because of Async copy
        amrex::ULong cpu_seed = 12345;
        amrex::ULong gpu_seed = 67890;
        amrex::InitRandom(cpu_seed, amrex::ParallelDescriptor::NProcs(),
                          gpu_seed);

        std::vector<SmallDataIO::column_t> write_data(3);
        for (auto &column : write_data)
        {
            column.resize(Npts);
            column = generate_random_numbers(Npts);
        }

        amrex::Gpu::streamSynchronize();

        static constexpr int data_precision = 10;
        static const amrex::Real err_tol = std::pow(10., -1.0 * data_precision);

        // Test the file writer: write out the random coordinates
        test_small_data_io_writer(write_data, data_precision);

        // Test the file reader: read the random numbers back in using their
        // column names

        // Columns do not have to be in order, this is tested below using
        // {"z", "x"}. But the first column returned will always go
        // with the first entry, the second column returned  with the second
        // entered etc.

        std::vector<std::string> read_these_columns{"z", "x"};
        auto read_data1 = test_small_data_io_reader(read_these_columns);

        // Test the file reader: read the data by specifying the column numbers
        const int min_col = 1;
        const int max_col = 3;
        auto read_data2   = test_small_data_io_reader(min_col, max_col);

        check_almost_equal(read_data1[0], write_data[2], err_tol);
        check_almost_equal(read_data1[1], write_data[0], err_tol);
        check_almost_equal(read_data2[0], write_data[0], err_tol);
        check_almost_equal(read_data2[1], write_data[1], err_tol);
        check_almost_equal(read_data2[2], write_data[2], err_tol);
    }

    amrex::Finalize();
}
