/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SMALLDATAIO_HPP_
#define SMALLDATAIO_HPP_

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <regex>
#include <string>
#include <vector>

//! A class for reading and writing small data to a file in ASCII format.
/*!
    A class for reading and writing small data, usually 0D, 1D or 2D, in ASCII
    format. For an example on how to use it, see the WeylExtraction class.
*/
class SmallDataIO
{
  public:
    using column_t = std::vector<double>;

    //! Choose between appending data to the same file, writing to a new file
    //! at each timestep or reading a file.
    enum Mode
    {
        APPEND, // data is APPENDed to the same file at each timestep
        NEW,    // data is written to a NEW file at each timestep
        READ    // read data
    };

    // A struct for information about the structure of a SmallDataIO file
    struct file_structure_t
    {
        int num_blocks{0}; // a block is separated by 2 blank lines
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

    // Maximum allowed file size in bytes when reading
    static constexpr int max_file_size = 1024 * 1024 * 1024;

  protected:
    std::string m_filename;
    double m_dt;
    double m_time;
    double m_restart_time;
    int m_step{};
    Mode m_mode;
    bool m_first_step; // this should be set to true if this is the first
                       // timestep
    static const std::string s_default_file_extension;
    static constexpr int s_default_data_precision = 10;
    int m_data_precision;
    int m_data_width;
    double m_data_epsilon; //!< the maximum data precision error
    static constexpr int s_default_coords_precision = 7;
    int m_coords_precision;
    int m_coords_width;
    double m_coords_epsilon; //!< the maximum coords precision error
    static constexpr int s_default_filename_steps_width = 6;

    std::fstream m_file;

    // These used for reading SmallDataIO files
    [[maybe_unused]] std::string m_file_contents;
    [[maybe_unused]] file_structure_t m_file_structure;
    [[maybe_unused]] bool m_structure_defined;

    // Reads the entire file
    [[maybe_unused]] void read_file();

    // Parses the file and determines its structure
    [[maybe_unused]] void determine_file_structure();

  public:
    //! The full constructor (reading and writing)
    SmallDataIO(const std::string &a_filename_prefix, double a_dt,
                double a_time, double a_restart_time, Mode a_mode,
                bool a_first_step, const file_structure_t *a_file_structure,
                const std::string &a_file_extension = s_default_file_extension,
                int a_data_precision                = s_default_data_precision,
                int a_coords_precision     = s_default_coords_precision,
                int a_filename_steps_width = s_default_filename_steps_width);

    //! Constructors for writing (opens file)
    SmallDataIO(const std::string &a_filename_prefix, double a_dt,
                double a_time, double a_restart_time, Mode a_mode,
                bool a_first_step,
                const std::string &a_file_extension = s_default_file_extension,
                int a_data_precision                = s_default_data_precision,
                int a_coords_precision     = s_default_coords_precision,
                int a_filename_steps_width = s_default_filename_steps_width);

    //! Old constructor which assumes SmallDataIO is called in
    //! specificPostTimeStep
    SmallDataIO(const std::string &a_filename_prefix, double a_dt,
                double a_time, double a_restart_time, Mode a_mode,
                const std::string &a_file_extension = s_default_file_extension,
                int a_data_precision                = s_default_data_precision,
                int a_coords_precision     = s_default_coords_precision,
                int a_filename_steps_width = s_default_filename_steps_width);

    //! Constructors for reading, when m_time, m_dt, m_restart_time are
    //! irrelevant
    SmallDataIO(const std::string &a_filename_prefix,
                const std::string &a_file_extension = s_default_file_extension,
                int a_data_precision                = s_default_data_precision,
                int a_coords_precision = s_default_coords_precision);

    // This version accepts an argument for the file structure as well (if
    // known)
    SmallDataIO(const std::string &a_filename_prefix,
                const file_structure_t *a_file_structure,
                const std::string &a_file_extension = s_default_file_extension,
                int a_data_precision                = s_default_data_precision,
                int a_coords_precision = s_default_coords_precision);

    //! Destructor (closes file)
    ~SmallDataIO();

    // disable default copy/move constructors and assignment operators
    SmallDataIO(const SmallDataIO &)            = delete;
    SmallDataIO &operator=(const SmallDataIO &) = delete;
    SmallDataIO(SmallDataIO &&)                 = delete;
    SmallDataIO &operator=(SmallDataIO &&)      = delete;

    // ------------ Writing Functions ------------

    //! Writes a header_line
    //! Use this for 0D or 1D data, where the first column is either the time
    //! or another coordinate whose name should be provided in
    //! a_pre_header_string.
    void write_header_line(const std::vector<std::string> &a_header_strings,
                           const std::string &a_pre_header_string = "time");

    //! Writes a header line
    //! Use this for 1D or 2D data when the first two or more columns are
    //! coordinates whose names should be provided in the vector of strings
    //! a_pre_header_strings
    void
    write_header_line(const std::vector<std::string> &a_header_strings,
                      const std::vector<std::string> &a_pre_header_strings);

    //! Writes a data line
    //! Use this for 0D or 1D data, where the first column is either the time or
    //! another coordinate.
    void write_data_line(const std::vector<double> &a_data,
                         const double a_coord);

    //! Writes a data line for a specific time.
    void write_time_data_line(const std::vector<double> &a_data);

    //! Writes a data line
    //! Use this for 1D or 2D data when the first two or more columns are
    //! coordinates.
    void write_data_line(const std::vector<double> &a_data,
                         const std::vector<double> &a_coords = {});

    //! This just adds a double line break to the file.
    void line_break();

    //! if restarting from an earlier checkpoint file, this function removes
    //! any time data that will be replaced.
    void remove_duplicate_time_data(const bool keep_m_time_data = false);

    // ------------ Reading Functions ------------

    //! Get the data associated to specific coordinates from the file
    //! Note only the first line with the given coordinates is obtained
    void get_specific_data_line(std::vector<double> &a_out_data,
                                const std::vector<double> &a_coords);

    //! Get the data associated to a specific coordinate (e.g. time) from the
    //! file
    void get_specific_data_line(std::vector<double> &a_out_data,
                                const double a_coord);

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

    // ------------ Other Functions --------------

    // Only rank designated as the IOProcessor reads. This is a helper function
    // to redistribute data amongst all ranks
    static void broadcast_data(std::vector<column_t> &data);

    //! returns the full filename of a file created in NEW mode at time=a_time
    //! with dt=a_dt
    static std::string get_new_filename(
        const std::string &a_file_prefix, double a_dt, double a_time,
        const std::string &a_file_extension = s_default_file_extension,
        int a_filename_steps_width          = s_default_filename_steps_width);

    //! returns m_data_epsilon
    double get_data_epsilon() const;

    //! returns the default data_epsilon
    static double get_default_data_epsilon();

    //! returns m_coords_epsilon
    double get_coords_epsilon() const;

    //! returns the default coords epsilon
    static double get_default_coords_epsilon();
};

#endif /* SMALLDATAIO_HPP_ */
