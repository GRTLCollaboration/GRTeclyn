/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SMALLDATAIOREADER_HPP
#define SMALLDATAIOREADER_HPP

#include <fstream>
#include <ios>
#include <sstream>
#include <string>
#include <vector>

#include "AMReX.H"
#include "AMReX_BLassert.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

// A class to read files written using SmallDataIO.

class SmallDataIOReader
{
  public:
    using column_t = std::vector<double>;
    // A struct for information about the structure of a SmallDataIO file
    struct file_structure_t
    {
        int num_blocks; // a block is separated by 2 blank lines
        std::vector<std::streamoff>
            block_starts; // position offsets from the beginning of the file

        std::vector<int> num_header_rows;  // the number of header rows in
                                           // each block
        std::vector<int> num_data_rows;    // the number of data rows in each
                                           // block
        std::vector<int> num_data_columns; // number of data columns in each
                                           // block
        void clear();
    };

  private:
    std::string m_filename;
    std::ifstream m_file;
    std::string m_file_contents;
    file_structure_t m_file_structure;
    bool m_structure_defined;
    int m_rank; // only the AMReX IO rank does the reading

    // Reads the entire file
    void read_file(const std::string &a_filename);

    // Parses the file and determines its structure
    void determine_file_structure();

  public:
    // Constructor
    SmallDataIOReader();

    // Destructor
    ~SmallDataIOReader();

    // Opens the file and sets m_filename. Note that this does not determine the
    // file structure
    void open(const std::string &a_filename);

    // Closes the file
    void close();

    // Set structure if known already (e.g. same as another file already
    // determined)
    void set_file_structure(const file_structure_t &a_file_structure);

    // File struture getter
    const file_structure_t &get_file_structure() const;

    // Print file structure
    void print_file_structure() const;

    // Get an interval of columns (inclusive) from a block
    void get_columns(std::vector<column_t> &out, int a_min_column,
                     int a_max_column, int a_block = 0);

    // Get columns based on their names in the header
    // (does not assume continguous column numbers)
    void get_columns(std::vector<column_t> &out,
                     const std::vector<std::string> &column_names,
                     const int a_block = 0);

    // Get all data columns from a block
    void get_all_data_columns(std::vector<column_t> &out, int a_block = 0);

    // Get a single column from a block
    void get_column(std::vector<column_t> &out, int a_column, int a_block = 0);

    // Get same data column from all blocks
    void get_data_column_from_all_blocks(std::vector<column_t> &out,
                                         int a_data_column);

    // Returns a vector of numeric values from a header row
    void get_data_from_header(std::vector<double> &out, int a_header_row_number,
                              int a_block) const;

    // Returns a vector of strings from a header row
    void get_header_strings(std::vector<std::string> &header,
                            int a_block = 0) const;

    // Utility function to skip the header rows to start reading where the data
    // is located
    void skip_ahead(std::istringstream &file_stream, int nlines_to_skip,
                    int a_block) const;

    // Only rank designated as the IOProcessor reads. This is a helper function
    // to redistribute data amongst all ranks
    static void broadcast_data(std::vector<SmallDataIOReader::column_t> &data);

    // Maximum allowed file size in bytes
    static constexpr int max_file_size = 1024 * 1024 * 1024;
};

#endif /* SMALLDATAIOREADER_HPP */
