/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "SmallDataIOReader.hpp"

#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <limits>
#include <regex>

void SmallDataIOReader::file_structure_t::clear()
{
    num_blocks = 0;
    block_starts.clear();
    num_data_rows.clear();
    num_header_rows.clear();
    num_data_columns.clear();
}

// Constructor
SmallDataIOReader::SmallDataIOReader()
    : m_structure_defined{false}, m_rank{amrex::ParallelDescriptor::MyProc()}
{
}

// Destructor
SmallDataIOReader::~SmallDataIOReader()
{
    close();
    if (m_file.is_open())
    {
        m_file.close();
    }
}

// Opens the file and sets m_filename and m_file_structure
void SmallDataIOReader::open(const std::string &a_filename)
{

    // TODO: move into constructor
    m_filename = a_filename;

    if (amrex::ParallelDescriptor::IOProcessor())
    {
        read_file(a_filename);
        determine_file_structure();
    }
}

// Closes the file
void SmallDataIOReader::close()
{
    m_filename.clear();
    m_structure_defined = false;

    if (amrex::ParallelDescriptor::IOProcessor())
    {
        m_file_structure.clear();
        m_file_contents.clear();
    }
}

void SmallDataIOReader::read_file(const std::string &a_filename)
{
    std::string warning{"Failed to open file: " + a_filename};

    if (std::ifstream input_stream{a_filename,
                                   std::ios::binary | std::ios::ate})
    {
        auto size = input_stream.tellg();

        AMREX_ASSERT_WITH_MESSAGE(
            size < max_file_size,
            "File size is too large for SmallDataIOReader.\n");

        std::string file_contents(size,
                                  '\0'); // construct string to stream size
        input_stream.seekg(0);
        input_stream.read(file_contents.data(), size);
        AMREX_ASSERT_WITH_MESSAGE(!file_contents.empty(),
                                  "File contents of " + a_filename +
                                      " could not be read");

        m_file_contents = file_contents;
    }
    else
    {
        amrex::Abort(warning);
    }
}

// Parses the file and determines its structure
// Don't call this directly, it will automatically be called when you open a
// file
void SmallDataIOReader::determine_file_structure()
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
void SmallDataIOReader::set_file_structure(
    const SmallDataIOReader::file_structure_t &a_file_structure)
{
    m_file_structure    = a_file_structure;
    m_structure_defined = true;
}

// File struture getter
const SmallDataIOReader::file_structure_t &
SmallDataIOReader::get_file_structure() const
{
    return m_file_structure;
}

// Utility for viewing the file struture

void SmallDataIOReader::print_file_structure() const
{
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        assert(m_structure_defined);
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
void SmallDataIOReader::get_columns(
    std::vector<SmallDataIOReader::column_t> &out, int a_min_column,
    int a_max_column, int a_block)
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

        for (int irow = 0; irow < m_file_structure.num_data_rows[a_block];
             ++irow)
        {
            double discard = 0.0;

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
void SmallDataIOReader::get_columns(
    std::vector<SmallDataIOReader::column_t> &out,
    const std::vector<std::string> &column_names, const int a_block)
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

        std::vector<int> column_numbers;

        for (const auto &column_name : column_names)
        {
            int column_no = 0;
            for (const auto &header_name : header)
            {
                if (header_name == column_name)
                {

                    column_numbers.push_back(column_no);
                    break;
                }
                column_no++;
            }

            // If all header_names have been checked but column_name has not
            // been found...
            if (column_no == m_file_structure.num_data_columns[a_block])
            {
                std::string error_message{
                    column_name + " could not be read. Please double check "
                                  "your inputs to SmallDataIOReader"};
                amrex::Abort(error_message);
            }
        }

        std::istringstream file_stream(m_file_contents);

        skip_ahead(file_stream, m_file_structure.num_header_rows[a_block],
                   a_block);

        for (int irow = 0; irow < m_file_structure.num_data_rows[a_block];
             ++irow)
        {
            double discard = 0.0;

            // This is a bit slow but returns the columns in the order that the
            // user requested them

            for (int icolumn = 0;
                 icolumn < m_file_structure.num_data_columns[a_block];
                 ++icolumn)
            {
                auto const column_ptr = std::find(
                    column_numbers.begin(), column_numbers.end(), icolumn);
                if (column_ptr != column_numbers.end())
                {
                    auto column_index = column_ptr - column_numbers.begin();
                    file_stream >> out[column_index][irow];
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
void SmallDataIOReader::get_all_data_columns(
    std::vector<SmallDataIOReader::column_t> &out, int a_block)
{
    assert(m_structure_defined);
    int min_data_column = 0;
    int max_data_column = m_file_structure.num_data_columns[a_block];
    get_columns(out, min_data_column, max_data_column, a_block);
}

void SmallDataIOReader::get_column(
    std::vector<SmallDataIOReader::column_t> &out, int a_column, int a_block)
{
    get_columns(out, a_column, a_column, a_block);
}

// Returns a vector of numeric values from a header row
// NOLINTBEGIN(bugprone-easily-swappable-parameters)
void SmallDataIOReader::get_data_from_header(std::vector<amrex::Real> &out,
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
void SmallDataIOReader::get_header_strings(std::vector<std::string> &header,
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

        //    std::vector<std::string>
        //    out(m_file_structure.num_data_columns[a_block]);

        for (auto &icol : header)
        {
            file_stream >> icol;
        }
    }
}
// Helper function to skip lines e.g. header lines
// NOLINTBEGIN(bugprone-easily-swappable-parameters)
void SmallDataIOReader::skip_ahead(std::istringstream &file_stream,
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

// Helper function to redistribute data amongst all ranks
void SmallDataIOReader::broadcast_data(
    std::vector<SmallDataIOReader::column_t> &data)
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

    amrex::ParallelDescriptor::Barrier();
}
