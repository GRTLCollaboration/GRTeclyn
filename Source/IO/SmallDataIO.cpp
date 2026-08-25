/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// Other includes
#include "SmallDataIO.hpp"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <AMReX_Utility.H>
#include <AMReX_Vector.H>

#include <cmath>
#include <random>
// (MR): if it were up to me, I'd be using the C++17 filesystems library
// instead of cstdio but I'm sure someone would tell me off for not maintaining
// backwards compatability.
#include <cstdio> // for std::rename and std::remove
#include <sstream>
#include <utility>

// ------------ Constructors -----------------

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
// This has to be initialised outside the class declaration in C++14
const std::string SmallDataIO::s_default_file_extension = ".dat";

SmallDataIO::SmallDataIO(const std::string &a_filename_prefix, amrex::Real a_dt,
                         amrex::Real a_time, amrex::Real a_restart_time,
                         Mode a_mode, bool a_first_step,
                         const file_structure_t *a_file_structure,
                         const std::string &a_file_extension,
                         int a_data_precision, int a_coords_precision,
                         int a_filename_steps_width)
    : m_filename(a_filename_prefix + a_file_extension), m_dt(a_dt),
      m_time(a_time), m_restart_time(a_restart_time), m_mode(a_mode),
      m_first_step(a_first_step), m_data_precision(a_data_precision),
      // data columns need extra space for scientific notation
      // compared to coords columns
      m_data_width(m_data_precision + 10),
      m_data_epsilon(std::pow(10.0, -a_data_precision)),
      m_coords_precision(a_coords_precision),
      m_coords_width(m_coords_precision + 5),
      m_coords_epsilon(std::pow(10.0, -a_coords_precision)),
      // this is for reading only
      m_structure_defined{false}
{

    if (amrex::ParallelDescriptor::IOProcessor())
    {
        std::ios::openmode file_openmode = std::ios::out;
        if (m_mode == APPEND)
        {
            if (m_first_step)
            {
                // overwrite any existing file if this is the first step
                file_openmode = std::ios::out;
            }
            else if (m_restart_time > 0. &&
                     m_time < m_restart_time + m_dt + m_coords_epsilon)
            {
                // allow reading in the restart case so that duplicate time
                // data may be removed
                file_openmode = std::ios::app | std::ios::in;
            }
            else
            {
                // default mode is just appending to existing file
                file_openmode = std::ios::app;
            }
        }
        else if (m_mode == NEW)
        {
            file_openmode = std::ios::out;
            m_filename =
                get_new_filename(a_filename_prefix, m_dt, m_time,
                                 a_file_extension, a_filename_steps_width);
        }
        else if (m_mode == READ)
        {
            file_openmode = std::ios::in;
        }
        else
        {
            amrex::Abort("SmallDataIO: mode not supported");
        }
        if (m_mode == APPEND && m_first_step || m_mode == NEW)
        {
            // Rather than overwriting files from previous simulations, we
            // rename the old files to "filename.old.<random string>" like AMReX
            // does for checkpoints and plotfiles
            bool call_mpi_barrier = false;
            // Even though "directory" is in this function name, it works fine
            // for any type of file.
            amrex::UtilRenameDirectoryToOld(m_filename, call_mpi_barrier);
        }
        m_file.open(m_filename, file_openmode);

        if (!m_file)
        {
            amrex::Abort("SmallDataIO: error opening file " + m_filename);
        }
        if (m_mode == READ)
        {
            read_file();
            if (a_file_structure == nullptr)
            {
                determine_file_structure();
            }
            else
            {
                set_file_structure(*a_file_structure);
            }
        }
    }
}

SmallDataIO::SmallDataIO(const std::string &a_filename_prefix, amrex::Real a_dt,
                         amrex::Real a_time, amrex::Real a_restart_time,
                         Mode a_mode, bool a_first_step,
                         const std::string &a_file_extension,
                         int a_data_precision, int a_coords_precision,
                         int a_filename_steps_width)
    : SmallDataIO(a_filename_prefix, a_dt, a_time, a_restart_time, a_mode,
                  a_first_step, nullptr, a_file_extension, a_data_precision,
                  a_coords_precision, a_filename_steps_width)
{
}

SmallDataIO::SmallDataIO(const std::string &a_filename_prefix, amrex::Real a_dt,
                         amrex::Real a_time, amrex::Real a_restart_time,
                         Mode a_mode, const std::string &a_file_extension,
                         int a_data_precision, int a_coords_precision,
                         int a_filename_steps_width)
    : SmallDataIO(a_filename_prefix, a_dt, a_time, a_restart_time, a_mode,
                  (a_time == a_dt), nullptr, a_file_extension, a_data_precision,
                  a_coords_precision, a_filename_steps_width)
{
}

SmallDataIO::SmallDataIO(const std::string &a_filename_prefix,
                         const std::string &a_file_extension,
                         int a_data_precision, int a_coords_precision)
    : SmallDataIO(a_filename_prefix, 0.0, 0.0, 0.0, READ, false, nullptr,
                  a_file_extension, a_data_precision, a_coords_precision, 0)
{
}

SmallDataIO::SmallDataIO(const std::string &a_filename_prefix,
                         const file_structure_t *a_file_structure,
                         const std::string &a_file_extension,
                         int a_data_precision, int a_coords_precision)
    : SmallDataIO(a_filename_prefix, 0.0, 0.0, 0.0, READ, false,
                  a_file_structure, a_file_extension, a_data_precision,
                  a_coords_precision, 0)
{
}

//! Destructor (closes file)
SmallDataIO::~SmallDataIO()
{
    m_filename.clear();
    m_structure_defined = false;

    if (amrex::ParallelDescriptor::IOProcessor())
    {
        m_file.close();
        m_file_structure.clear();
        m_file_contents.clear();
    }
}

// NOLINTEND(bugprone-easily-swappable-parameters)

// ------------ Writing Functions ------------

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
void SmallDataIO::write_header_line(
    const std::vector<std::string> &a_header_strings,
    const std::string &a_pre_header_string)
{
    std::vector<std::string> pre_header_strings;
    if (!a_pre_header_string.empty())
    {
        pre_header_strings.push_back(a_pre_header_string);
    }
    write_header_line(a_header_strings, pre_header_strings);
}

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
void SmallDataIO::write_header_line(
    const std::vector<std::string> &a_header_strings,
    const std::vector<std::string> &a_pre_header_strings)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        // all header lines start with a '#'.
        m_file << "#";
        for (std::size_t istr = 0; istr < a_pre_header_strings.size(); ++istr)
        {
            // first column header is shorter due to preceeding #
            if (istr == 0)
            {
                m_file << std::setw(m_coords_width - 1)
                       << a_pre_header_strings[istr];
            }
            else
            {
                m_file << std::setw(m_coords_width)
                       << a_pre_header_strings[istr];
            }
        }
        for (const std::string &header_item : a_header_strings)
        {
            m_file << std::setw(m_data_width) << header_item;
        }
        m_file << "\n";
    }
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
void SmallDataIO::write_data_line(const std::vector<amrex::Real> &a_data,
                                  const amrex::Real a_coord)
{
    const std::vector<amrex::Real> coords(1, a_coord);
    write_data_line(a_data, coords);
}

void SmallDataIO::write_time_data_line(const std::vector<amrex::Real> &a_data)
{
    write_data_line(a_data, m_time);
}

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
void SmallDataIO::write_data_line(const std::vector<amrex::Real> &a_data,
                                  const std::vector<amrex::Real> &a_coords)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        m_file << std::fixed << std::setprecision(m_coords_precision);
        for (amrex::Real coord : a_coords)
        {
            m_file << std::setw(m_coords_width) << coord;
        }
        m_file << std::scientific << std::setprecision(m_data_precision);
        for (amrex::Real data : a_data)
        {
            m_file << std::setw(m_data_width) << data;
        }
        m_file << "\n";
    }
}

void SmallDataIO::line_break()
{
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        m_file << "\n\n";
    }
}

void SmallDataIO::remove_duplicate_time_data(const bool keep_m_time_data)
{
    if (amrex::ParallelDescriptor::IOProcessor() && m_restart_time > 0. &&
        m_mode == APPEND && m_time < m_restart_time + m_dt + m_coords_epsilon)
    {
        // copy lines with time < m_time into a temporary file
        m_file.seekg(0);
        std::string line;
        // adding a random integer might make this a little more robust...
        const unsigned long random_int = std::default_random_engine()();
        std::string temp_filename =
            m_filename + ".temp" + std::to_string(random_int);
        std::ofstream temp_file(temp_filename);
        int sign = -1;
        if (keep_m_time_data)
        {
            sign = 1;
        }
        while (std::getline(m_file, line))
        {
            if (!(line.find('#') == std::string::npos) ||
                std::stod(line.substr(0, m_coords_width)) <
                    m_time + sign * m_coords_epsilon)
            {
                temp_file << line << "\n";
            }
        }

        m_file.close();
        temp_file.close();

        // now delete the original file and rename the temporary file with the
        // original filename
        std::remove(m_filename.data());
        std::rename(temp_filename.data(), m_filename.data());
        // reopen the file in append mode
        m_file.open(m_filename, std::ios::app);
    }
}

// ------------ Reading Functions ------------

void SmallDataIO::file_structure_t::clear()
{
    num_blocks = 0;
    block_starts.clear();
    num_data_rows.clear();
    num_header_rows.clear();
    num_data_columns.clear();
}

// Only for use in the constructor - this assumes that the file is already open
void SmallDataIO::read_file()
{

    AMREX_ASSERT_WITH_MESSAGE(
        m_file, "This function can only be called from the constructor! \n");

    m_file.seekg(0, std::ios::end);

    auto size = m_file.tellg();

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        size < max_file_size,
        "SmallDataIO: File size is too large. Max file size is 1 GB. \n");

    std::string file_contents(size,
                              '\0'); // construct string to stream size
    m_file.seekg(0);
    m_file.read(file_contents.data(), size);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!file_contents.empty(),
                                     "SmallDataIO: File contents of " +
                                         m_filename + " could not be read");

    m_file_contents = file_contents;
}

// Parses the file and determines its structure
// Don't call this directly, it will automatically be called when you open a
// file
void SmallDataIO::determine_file_structure()
{
    // first check that the file has been read
    AMREX_ASSERT(m_file_contents.empty() == false);

    std::istringstream file_stream(m_file_contents);
    // go through each line and determine structure
    std::string line;
    int block_counter        = 0; // assume we always have one block
    auto current_position    = file_stream.tellg();
    int block_start_position = current_position;
    int header_row_counter   = 0;
    int data_row_counter     = 0;

    while (std::getline(file_stream, line))
    {
        if (!line.empty())
        {
            // header rows start with '#'
            if (line.find('#') == 0)
            {
                if (header_row_counter == 0)
                {
                    block_start_position = current_position;
                }

                ++header_row_counter;
            }
            else
            {
                if (data_row_counter++ == 0)
                {

                    // only count a new block if it contains a data row
                    m_file_structure.block_starts.push_back(
                        block_start_position);
                    ++block_counter;
                    // determine column structure from first data row in
                    // block get a vector of the widths of the columns
                    // including preceeding whitespace
                    std::string::size_type start_whitespace = 0;

                    int ncols{0};
                    while (!(start_whitespace == std::string::npos))
                    {
                        std::string::size_type start_non_whitespace =
                            line.find_first_not_of(' ', start_whitespace);
                        std::string::size_type next_start_whitespace =
                            line.find_first_of(' ', start_non_whitespace);
                        start_whitespace = next_start_whitespace;
                        ncols++;
                    }

                    m_file_structure.num_data_columns.push_back(ncols);
                }
            }
        }
        else
        {
            if (header_row_counter > 0 || data_row_counter > 0)

            {
                // end of previous block
                m_file_structure.num_header_rows.push_back(header_row_counter);
                m_file_structure.num_data_rows.push_back(data_row_counter);
                // reset the counters
                header_row_counter = 0;
                data_row_counter   = 0;
            }
        }
        current_position = file_stream.tellg();
    }
    // Just in case the file ends without a line break:
    if (data_row_counter > 0)
    {
        m_file_structure.num_header_rows.push_back(header_row_counter);
        m_file_structure.num_data_rows.push_back(data_row_counter);
    }

    m_file_structure.num_blocks = block_counter;

    AMREX_ASSERT(m_file_structure.num_data_rows.size() ==
                 m_file_structure.num_blocks);
    AMREX_ASSERT(m_file_structure.num_header_rows.size() ==
                 m_file_structure.num_blocks);
    AMREX_ASSERT(m_file_structure.num_data_columns.size() ==
                 m_file_structure.num_blocks);
    AMREX_ASSERT(m_file_structure.block_starts.size() ==
                 m_file_structure.num_blocks);

    m_structure_defined = true;
}

// Set structure if known already (e.g. same as another file already
// determined)
void SmallDataIO::set_file_structure(
    const SmallDataIO::file_structure_t &a_file_structure)
{
    m_file_structure    = a_file_structure;
    m_structure_defined = true;
}

// File struture getter
const SmallDataIO::file_structure_t &SmallDataIO::get_file_structure() const
{
    return m_file_structure;
}

// Utility for viewing the file struture

void SmallDataIO::print_file_structure() const
{
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        AMREX_ASSERT(m_structure_defined);

        amrex::Print() << "Total number of blocks: "
                       << m_file_structure.num_blocks << "\n";

        for (int i = 0; i < m_file_structure.num_blocks; i++)
        {

            // Printing info on file structure:
            amrex::Print() << "#######" << "\n";
            amrex::Print() << "Info on Block " << i << ":" << "\n";
            amrex::Print() << "Block start: "
                           << m_file_structure.block_starts[i] << "\n";
            amrex::Print() << "Number of columns in Block " << i << ": "
                           << m_file_structure.num_data_columns[i] << "\n";
            amrex::Print() << "Number of rows in Block " << i << ": "
                           << m_file_structure.num_data_rows[i] << "\n";
            amrex::Print() << "Number of header rows in Block " << i << ": "
                           << m_file_structure.num_header_rows[i] << "\n";

            // Printing column names:

            std::vector<std::string> header;
            get_header_strings(header, i);

            amrex::Print() << "Header strings found: ";
            for (auto &column_name : header)
            {
                amrex::Print() << column_name << " ";
            }
            amrex::Print() << "\n";
            amrex::Print() << "#######" << "\n";
        }
    }
}

// Get an interval of columns (inclusive) from a block
// NOLINTBEGIN(bugprone-easily-swappable-parameters)
void SmallDataIO::get_columns(std::vector<SmallDataIO::column_t> &out,
                              int a_min_column, int a_max_column, int a_block)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        AMREX_ASSERT(m_file_contents.empty() == false);
        AMREX_ASSERT(m_structure_defined);
        AMREX_ASSERT(0 <= a_min_column);
        AMREX_ASSERT(a_max_column <=
                     m_file_structure.num_data_columns[a_block]);
        AMREX_ASSERT(a_min_column <= a_max_column);
        const int num_columns = a_max_column - a_min_column + 1;

        out.resize(num_columns);
        for (auto &column : out)
        {
            column.resize(m_file_structure.num_data_rows[a_block]);
        }

        std::istringstream file_stream(m_file_contents);
        skip_ahead(file_stream, m_file_structure.num_header_rows[a_block],
                   a_block);

        amrex::Real discard = 0.0;
        for (int irow = 0; irow < m_file_structure.num_data_rows[a_block];
             ++irow)
        {

            for (int icolumn = 0;
                 icolumn < m_file_structure.num_data_columns[a_block];
                 ++icolumn)
            {
                if (a_min_column <= icolumn && icolumn <= a_max_column)
                {

                    file_stream >> out[icolumn - a_min_column][irow];
                }
                else
                {
                    file_stream >> discard;
                }
            }
        }
    }
}

// Get columns from a block based on their names in the header
// NOLINTBEGIN(bugprone-easily-swappable-parameters)
void SmallDataIO::get_columns(std::vector<SmallDataIO::column_t> &out,
                              const std::vector<std::string> &column_names,
                              const int a_block)
// NOLINTEND(bugprone-easily-swappable-parameters)
{

    if (amrex::ParallelDescriptor::IOProcessor())
    {
        // Resize for the number of requested columns
        out.resize(column_names.size());

        // Resize for the total number of rows
        for (auto &column : out)
        {
            column.resize(m_file_structure.num_data_rows[a_block]);
        }

        AMREX_ASSERT(m_file_contents.empty() == false);
        AMREX_ASSERT(m_structure_defined);

        // Get column names from header
        // Assume that the column names are in the line immediately proceeding
        // the data, like this:
        // # blah blah this is my data blah blah
        // # blah blah more things to say blah blah
        // # x y z

        std::vector<std::string> header;
        get_header_strings(header, a_block);

        // Associate each column name with the number
        // of where they are in the file

        // First entry in the std::map is the column number in the file
        // (file_column) Second entry in the std::map is the column number given
        // by the user (input_column)
        std::map<int, int> column_mapping;

        bool found{false};

        for (int input_column = 0; input_column < column_names.size();
             input_column++)
        {
            for (int file_column = 0; file_column < header.size();
                 file_column++)
            {
                if (header[file_column] == column_names[input_column])
                {
                    column_mapping.emplace(file_column, input_column);
                    found = true;
                }
            }
            // If all header_names have been checked but column_name has not
            // been found...
            if (!found)
            {
                std::string error_message{
                    "SmallDataIO: " + column_names[input_column] +
                    " could not be read. Please double check "
                    "your inputs to SmallDataIO::get_columns"};
                amrex::Abort(error_message);
            }
        }

        const int ncols = m_file_structure.num_data_columns[a_block];

        // To avoid calling std::map::find for each row of data,
        // we define a bool for each column (in the original file ordering) so
        // we know if that value needs to be saved or not.
        std::vector<bool> columns_requested(ncols, false);

        for (auto const &[file_order, user_order] : column_mapping)
        {
            columns_requested[file_order] = true;
        }

        std::istringstream file_stream(m_file_contents);

        skip_ahead(file_stream, m_file_structure.num_header_rows[a_block],
                   a_block);

        amrex::Real discard = 0.0;
        for (int row = 0; row < m_file_structure.num_data_rows[a_block]; ++row)
        {
            for (int file_column = 0;
                 file_column < m_file_structure.num_data_columns[a_block];
                 ++file_column)
            {
                if (columns_requested[file_column])
                {
                    // The mapping is there to return the columns in the same
                    // order as the user input

                    file_stream >> out[column_mapping.at(file_column)][row];
                }
                else
                {
                    file_stream >> discard;
                }
            }
        }
    }
}

// Get a data column from a block
void SmallDataIO::get_all_data_columns(std::vector<SmallDataIO::column_t> &out,
                                       int a_block)
{
    assert(m_structure_defined);
    int min_data_column = 0;
    int max_data_column = m_file_structure.num_data_columns[a_block];
    get_columns(out, min_data_column, max_data_column, a_block);
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
void SmallDataIO::get_column(std::vector<SmallDataIO::column_t> &out,
                             int a_column, int a_block)
{
    get_columns(out, a_column, a_column, a_block);
}

// Returns a vector of numeric values from a header row
// NOLINTBEGIN(bugprone-easily-swappable-parameters)
void SmallDataIO::get_data_from_header(std::vector<amrex::Real> &out,
                                       int a_header_row_number,
                                       int a_block) const
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    if (amrex::ParallelDescriptor::IOProcessor())
    {

        assert(m_file_contents.empty() == false);
        assert(m_structure_defined);
        assert(a_header_row_number < m_file_structure.num_header_rows[a_block]);

        std::istringstream file_stream(m_file_contents);
        skip_ahead(file_stream, a_header_row_number, a_block);

        // get desired header line
        std::string line;
        std::getline(file_stream, line);

        // find numbers in header using regex
        // I think this takes a long time to compile...
        // NOLINTNEXTLINE(bugprone-unused-local-non-trivial-variable)
        std::regex number("[+-]?([0-9]*\\.)?[0-9]+");
        auto numbers_begin =
            std::sregex_iterator(line.begin(), line.end(), number);
        auto numbers_end = std::sregex_iterator();

        for (std::sregex_iterator rit = numbers_begin; rit != numbers_end;
             ++rit)
        {
            // put matches in vector
            out.push_back(std::stod((*rit).str()));
        }
    }
}

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
void SmallDataIO::get_header_strings(std::vector<std::string> &header,
                                     int a_block) const
// NOLINTEND(bugprone-easily-swappable-parameters)
{

    if (amrex::ParallelDescriptor::IOProcessor())
    {
        assert(m_file_contents.empty() == false);
        assert(m_structure_defined);

        header.resize(m_file_structure.num_data_columns[0]);

        std::istringstream file_stream(m_file_contents);
        int nlines_to_skip = m_file_structure.num_header_rows[a_block] - 1;
        skip_ahead(file_stream, nlines_to_skip, a_block);

        std::string hashes; // remove the #s at the beginning of the line
        file_stream >> hashes;

        for (auto &icol : header)
        {
            file_stream >> icol;
        }
    }
}
// Helper function to skip lines e.g. header lines
// NOLINTBEGIN(bugprone-easily-swappable-parameters)
void SmallDataIO::skip_ahead(std::istringstream &file_stream,
                             int nlines_to_skip, int a_block) const
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    // rewind to the beginning of the block
    file_stream.seekg(m_file_structure.block_starts[a_block], std::ios::beg);

    // assume header rows are at start of block
    for (int irow = 0; irow < nlines_to_skip; ++irow)
    {
        // skip lines before desired row
        file_stream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
}

void SmallDataIO::get_specific_data_line(
    std::vector<amrex::Real> &a_out_data,
    const std::vector<amrex::Real> &a_coords)
{
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        bool line_found = false;
        // first set the current position to the beginning of the file
        m_file.seekg(0);

        // get a string of the coords as they are in the file
        std::stringstream coords_ss;
        coords_ss << std::fixed << std::setprecision(m_coords_precision);
        for (amrex::Real coord : a_coords)
        {
            coords_ss << std::setw(m_coords_width) << coord;
        }
        std::string coords_string = coords_ss.str();

        // now search for lines that start with coords_string and put the data
        // in a_out_data
        std::string line;
        while (std::getline(m_file, line))
        {
            if (!(line.find(coords_string) == std::string::npos))
            {
                for (std::size_t ichar = a_coords.size() * m_coords_width;
                     ichar < line.size(); ichar += m_data_width)
                {
                    amrex::Real data_value =
                        std::stod(line.substr(ichar, m_data_width));
                    a_out_data.push_back(data_value);
                }
                line_found = true;
                // only want the first occurrence so break the while loop
                break;
            }
        }
        if (!line_found)
        {
            amrex::Abort(
                "SmallDataIO : Data to be read in at coord not found in file");
        }
    }
    // Optionally, broadcast the data using SmallDataIO::broadcast_data
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
void SmallDataIO::get_specific_data_line(std::vector<amrex::Real> &a_out_data,
                                         const amrex::Real a_coord)
{
    std::vector<amrex::Real> coords(1, a_coord);
    get_specific_data_line(a_out_data, coords);
}

// ------------ Other Functions --------------

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
std::string SmallDataIO::get_new_filename(const std::string &a_file_prefix,
                                          amrex::Real a_dt, amrex::Real a_time,
                                          const std::string &a_file_extension,
                                          int a_filename_steps_width)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    AMREX_ASSERT(a_dt > 0);
    const int step = static_cast<int>(std::round(a_time / a_dt));

    // append step number to filename (pad to make it
    // a_filename_steps_width digits).
    std::string step_string = std::to_string(step);
    if (a_filename_steps_width < static_cast<int>(step_string.length()))
    {
        amrex::Abort("SmallDataIO: a_filename_steps_width too small "
                     "for step number");
    }
    std::string step_string_padded =
        std::string(a_filename_steps_width - step_string.length(), '0') +
        step_string;
    // append step number to filename if in NEW mode
    return a_file_prefix + step_string_padded + a_file_extension;
}

// returns m_data_epsilon
amrex::Real SmallDataIO::get_data_epsilon() const { return m_data_epsilon; }

// returns the default data_epsilon
amrex::Real SmallDataIO::get_default_data_epsilon()
{
    return pow(10.0, -s_default_data_precision);
}

// returns m_coords_epsilon
amrex::Real SmallDataIO::get_coords_epsilon() const { return m_coords_epsilon; }

// returns the default coords epsilon
amrex::Real SmallDataIO::get_default_coords_epsilon()
{
    return pow(10.0, -s_default_coords_precision);
}

// Helper function to redistribute data amongst all ranks
// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
void SmallDataIO::broadcast_data(std::vector<SmallDataIO::column_t> &data)
{
    int nrows{0};
    int ncols{0};

    if (amrex::ParallelDescriptor::IOProcessor())
    {
        nrows = data[0].size();
        ncols = data.size();
    }
    amrex::ParallelDescriptor::Bcast(
        &ncols, 1, amrex::ParallelDescriptor::IOProcessorNumber());

    amrex::ParallelDescriptor::Bcast(
        &nrows, 1, amrex::ParallelDescriptor::IOProcessorNumber());

    data.resize(ncols);

    for (auto &column : data)
    {
        column.resize(nrows);
        amrex::ParallelDescriptor::Bcast(
            column.data(), nrows,
            amrex::ParallelDescriptor::IOProcessorNumber());
    }
}
