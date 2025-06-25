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
SmallDataIOReader::SmallDataIOReader() : m_structure_defined{false} {}

// Destructor
SmallDataIOReader::~SmallDataIOReader()
{
    if (m_file.is_open())
    {
        close();
    }
}

// Opens the file and sets m_filename. Note that this does not determine the
// file structure
void SmallDataIOReader::open(const std::string &a_filename)
{

    // TODO: move into constructor

    m_filename          = a_filename;
    m_structure_defined = false;
    m_file_contents     = read_entire_file(a_filename);

    //    m_file.open(m_filename);

    // check file opening successful
    // if (!m_file)
    // {
    //     std::cerr << "Error in opening " << m_filename << ". Exiting..."
    //               << std::endl;
    //     exit(1);
    // }
}

// Closes the file
void SmallDataIOReader::close()
{
    // TODO: move into destructor

    // if (m_file.is_open())
    // {
    //     m_file.close();
    // }
    m_filename.clear();
    m_structure_defined = false;
    m_file_structure.clear();
    m_file_contents.clear();
}

std::string SmallDataIOReader::read_entire_file(const std::string &a_filename)
{

    std::string warning{"Failed to open file \n"};

    if (std::ifstream input_stream{a_filename,
                                   std::ios::binary | std::ios::ate})
    {
        auto size = input_stream.tellg();
        if (size > max_file_size)
        {
            return ("File size is too large for SmallDataIOReader.\n");
            exit(1);
        }

        std::string file_contents(size,
                                  '\0'); // construct string to stream size
        input_stream.seekg(0);
        if (input_stream.read(file_contents.data(), size))
        {
            //            std::cout << file_contents << '\n';
            return (file_contents);
        }
        else
        {
            return (warning);
            exit(1);
        }
    }
    else
    {
        return (warning);
        exit(1);
    }
}

// Parses the file and determines its structure
void SmallDataIOReader::determine_file_structure()
{
    // first check that the file has been read
    assert(m_file_contents.empty() == false);

    // move file stream position to start of file
    // m_file.clear();
    // m_file.seekg(0, std::ios::beg);

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
                    // determine column structure from first data row in block
                    // get a vector of the widths of the columns including
                    // preceeding whitespace
                    std::string::size_type start_whitespace = 0;
                    std::vector<int> data_widths{};
                    while (!(start_whitespace == std::string::npos))
                    {
                        std::string::size_type start_non_whitespace =
                            line.find_first_not_of(' ', start_whitespace);
                        std::string::size_type next_start_whitespace =
                            line.find_first_of(' ', start_non_whitespace);
                        std::string::size_type width = 0;
                        if (next_start_whitespace == std::string::npos)
                        {
                            width = line.length() - start_whitespace;
                        }
                        else
                        {
                            width = next_start_whitespace - start_whitespace;
                        }
                        data_widths.push_back(width);
                        start_whitespace = next_start_whitespace;
                    }

                    // if (block_counter == 1)
                    // {
                    //     // first data row in file so get coord and data width
                    //     // from this. assume min width is coord width and max
                    //     is
                    //     // data width
                    //     auto widths_minmax_it = std::minmax_element(
                    //         m_file_structure.data_widths.begin(),
                    //         m_file_structure.data_widths.end());
                    //     m_file_structure.coords_width =
                    //         *(widths_minmax_it.first);
                    //     m_file_structure.data_width =
                    //         *(widths_minmax_it.second);
                    // }
                    m_file_structure.num_data_columns.push_back(
                        data_widths.size());
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
        header_row_counter = 0;
        data_row_counter   = 0;
    }

    m_file_structure.num_blocks = block_counter;

    assert(m_file_structure.num_data_rows.size() ==
           m_file_structure.num_blocks);
    assert(m_file_structure.num_header_rows.size() ==
           m_file_structure.num_blocks);
    assert(m_file_structure.num_data_columns.size() ==
           m_file_structure.num_blocks);
    assert(m_file_structure.block_starts.size() == m_file_structure.num_blocks);

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

// Get an interval of columns (inclusive) from a block
// NOLINTBEGIN(bugprone-easily-swappable-parameters)
std::vector<SmallDataIOReader::column_t>
SmallDataIOReader::get_columns(int a_min_column, int a_max_column, int a_block)
// NOLINTEND(bugprone-easily-swappable-parameters)
{

    assert(m_file_contents.empty() == false);
    assert(m_structure_defined);
    assert(0 <= a_min_column);
    assert(a_max_column <= m_file_structure.num_data_columns[a_block]);
    assert(a_min_column <= a_max_column);
    const int num_columns = a_max_column - a_min_column + 1;
    std::vector<column_t> out(num_columns);
    for (auto &column : out)
    {
        column.resize(m_file_structure.num_data_rows[a_block]);
    }

    // move stream position to start of block
    //    m_file_stream.clear();

    std::istringstream file_stream(m_file_contents);
    file_stream.seekg(m_file_structure.block_starts[a_block], std::ios::beg);

    // assume header rows are all at the top of the block so skip these
    for (int irow = 0; irow < m_file_structure.num_header_rows[a_block]; ++irow)
    {
        file_stream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    for (int irow = 0; irow < m_file_structure.num_data_rows[a_block]; ++irow)
    {
        //        std::getline(file_stream, line);
        //        std::istringstream ss(line);
        // std::string temp;
        // char del = ' ';
        double discard = 0.0;

        for (int icolumn = 0;
             icolumn < m_file_structure.num_data_columns[a_block]; ++icolumn)
        {
            if (a_min_column <= icolumn && icolumn < a_max_column)
            {

                file_stream >> out[icolumn - a_min_column][irow];
                std::cout << out[icolumn - a_min_column][irow] << '\n';
            }
            else
            {
                file_stream >> discard;
            }
        }
        std::cout << '\n';
    }

    return out;
}

// Get an interval of columns (inclusive) from a block
// NOLINTBEGIN(bugprone-easily-swappable-parameters)
std::vector<SmallDataIOReader::column_t>
SmallDataIOReader::get_columns(const std::string &column_names,
                               const int a_block)
// NOLINTEND(bugprone-easily-swappable-parameters)
{

    // in case there are several entries, split column_name into separate names
    std::istringstream name_stream(column_names);
    std::vector<std::string> separate_names;

    while (!name_stream.eof())
    {
        std::string name;
        name_stream >> name;
        separate_names.push_back(name);
    }

    // get column names from header
    auto header = get_header_strings(1, a_block);

    std::string test_string;
    for (const auto &piece : header)
    {
        test_string += piece;
        test_string += " ";
    }
    std::cout << test_string << '\n';

    std::vector<int> column_numbers;

    for (const auto &name : separate_names)
    {
        for (int j = 0; j < m_file_structure.num_data_columns[a_block]; ++j)
        {
            const std::string &current_name = header[j];
            if (current_name.find(name) == 0)
            {
                column_numbers.push_back(j);
                std::cout << "Found column name " << name << '\n';
                break;
            }
        }
        // auto col_num = test_string.find(indiv_names[i]);
        // std::cout << col_num << std::endl;
        // column_numbers.push_back(col_num);
    }
    // for (int i=0; i < m_file_structure.num_data_columns[a_block]; ++i)
    // {
    //   std::string header_name = header[i];

    //   //        std::cout << *it << std::endl;
    // 	if(header_name.find(indiv_names[0])==0)
    // 	  {
    // 	    column_numbers.push_back(i);
    // 	    std::cout << "Found column name" << std::endl;
    // 	    break;
    // 	  }
    // }

    std::sort(column_numbers.begin(), column_numbers.end());

    // for (auto it = column_numbers.begin(); it != column_numbers.end(); ++it)
    // {
    //     std::cout << *it << '\n';
    // }

    // placeholders for now
    auto want_cols =
        std::minmax_element(column_numbers.begin(), column_numbers.end());
    int a_min_column = *want_cols.first;
    int a_max_column = *want_cols.second + 1;

    assert(m_file_contents.empty() == false);
    assert(m_structure_defined);
    assert(0 <= a_min_column);
    assert(a_max_column <= m_file_structure.num_data_columns[a_block]);
    assert(a_min_column <= a_max_column);
    const int num_columns = a_max_column - a_min_column + 1;
    std::vector<column_t> out(num_columns);
    for (auto &column : out)
    {
        column.resize(m_file_structure.num_data_rows[a_block]);
    }

    // move stream position to start of block
    //    m_file_stream.clear();

    std::istringstream file_stream(m_file_contents);
    file_stream.seekg(m_file_structure.block_starts[a_block], std::ios::beg);

    // assume header rows are all at the top of the block so skip these
    for (int irow = 0; irow < m_file_structure.num_header_rows[a_block]; ++irow)
    {
        file_stream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    for (int irow = 0; irow < m_file_structure.num_data_rows[a_block]; ++irow)
    {
        //        std::getline(file_stream, line);
        //        std::istringstream ss(line);
        // std::string temp;
        // char del = ' ';
        double discard = 0.0;

        // skip over to minimum column number in this line
        for (int icolumn = 0; icolumn < a_min_column; ++icolumn)
        {
            file_stream >> discard;
        }

        int count = 0;
        for (int icolumn = a_min_column; icolumn < a_max_column; ++icolumn)
        {
            if (column_numbers[count] == icolumn)

            {

                file_stream >> out[count][irow];
                count++;
            }
            else
            {
                file_stream >> discard;
            }
        }

        // skip over max column number in this line
        for (int icolumn = a_max_column;
             icolumn < m_file_structure.num_data_columns[a_block]; ++icolumn)
        {
            file_stream >> discard;
        }

        // std::cout << std::endl;
    }

    return out;
}

// Get a data column from a block
std::vector<SmallDataIOReader::column_t>
SmallDataIOReader::get_all_data_columns(int a_block)
{
    assert(m_structure_defined);
    int min_data_column = 0;
    int max_data_column = m_file_structure.num_data_columns[a_block];
    return get_columns(min_data_column, max_data_column, a_block);
}

SmallDataIOReader::column_t SmallDataIOReader::get_column(int a_column,
                                                          int a_block)
{
    auto out_vect = get_columns(a_column, a_column, a_block);
    return out_vect[0];
}

// Returns a vector of numeric values from a header row
// NOLINTBEGIN(bugprone-easily-swappable-parameters)
std::vector<double>
SmallDataIOReader::get_data_from_header(int a_header_row_number, int a_block)
// NOLINTEND(bugprone-easily-swappable-parameters)
{

    assert(m_file_contents.empty() == false);
    assert(m_structure_defined);
    assert(a_header_row_number < m_file_structure.num_header_rows[a_block]);

    // move stream to start of block
    // m_file.clear();
    std::istringstream file_stream(m_file_contents);
    file_stream.seekg(m_file_structure.block_starts[a_block], std::ios::beg);
    std::string line;

    // assume header rows are at start of block
    for (int irow = 0; irow < a_header_row_number; ++irow)
    {
        // skip lines before desired row
        file_stream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    // get desired header line
    std::getline(file_stream, line);

    // find numbers in header using regex
    // I think this takes a long time to compile...
    std::regex number("[+-]?([0-9]*\\.)?[0-9]+");
    auto numbers_begin = std::sregex_iterator(line.begin(), line.end(), number);
    auto numbers_end   = std::sregex_iterator();
    std::vector<double> out;
    for (std::sregex_iterator rit = numbers_begin; rit != numbers_end; ++rit)
    {
        // put matches in vector
        out.push_back(std::stod((*rit).str()));
    }

    return out;
}
// NOLINTBEGIN(bugprone-easily-swappable-parameters)
std::vector<std::string>
SmallDataIOReader::get_header_strings(int a_header_row_number, int a_block)
// NOLINTEND(bugprone-easily-swappable-parameters)
{

    assert(m_file_contents.empty() == false);
    assert(m_structure_defined);
    assert(a_header_row_number < m_file_structure.num_header_rows[a_block]);

    // move stream to start of block
    //    m_file.clear();
    std::istringstream file_stream(m_file_contents);
    file_stream.seekg(m_file_structure.block_starts[a_block], std::ios::beg);

    // assume header rows are at start of block
    for (int irow = 0; irow < a_header_row_number; ++irow)
    {
        // skip lines before desired row
        file_stream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    std::string hashes; // remove the #s at the beginning of the line
    file_stream >> hashes;

    std::vector<std::string> out(m_file_structure.num_data_columns[a_block]);

    std::map<std::string, int> m;
    for (auto &icol : out)
    {
        file_stream >> icol;
        //	m.emplace(out[icol], icol);
    }

    // For (int icol = 0; icol < out.size(); ++icol)
    // {
    //     // int start_idx = m_file_structure.num_coords_columns[a_block] *
    //     //                     m_file_structure.coords_width +
    //     //                 icol * m_file_structure.data_width;
    //   int start_idx = icol * m_file_structure.data_width;
    //     if (start_idx < line.size())
    //     {
    //         std::string column_header =
    //             line.substr(start_idx, m_file_structure.data_width);
    //         int first_non_whitespace_char = std::distance(
    //             column_header.begin(),
    //             std::find_if(column_header.begin(), column_header.end(),
    //                          [](char c) { return (c != ' '); }));
    //         out[icol] = column_header.substr(first_non_whitespace_char);
    //     }
    //     else
    //     {
    //         out[icol] = "";
    //     }
    // }

    return out;
}
