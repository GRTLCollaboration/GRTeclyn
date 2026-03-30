/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef INFLATIONEXTRACTION_HPP_
#define INFLATIONEXTRACTION_HPP_

#include "InflationConfig.hpp"
#include "VarsTools.hpp"
#include "FilesystemTools.hpp" 
#include "SmallDataIO.hpp"

#include <AMReX_Print.H>
#include <AMReX_Vector.H>
#include <AMReX_FFT.H>

using namespace amrex;

class InflationExtraction
{
    public:
        // Names of diagnostic variables
        static inline const Vector<std::string> var_names = {"R", "hplus", "hcross"};

        // Constructor used in extraction of diagnostics
        InflationExtraction(InflationConfig a_params)
                : m_params(a_params), norm(std::pow(a_params.L, -3.))   
        {}

        // Main routines
        void derive(const MultiFab &source, MultiFab &out, const int dcomp);
        void extract(const MultiFab &state);

        // Set-up function when printing to data files
        void set_print_params(const std::string data_path, const Real dt,  
                     const Real cur_time, const int restart_time, const int first_step)
        {
            m_data_path = data_path;
            m_dt = dt;
            m_cur_time = cur_time;
            m_restart_time = restart_time;
            m_first_step = first_step;
        }

        Vector<Real> print_moment(MultiFab &field, const Vector<std::string> names,  
                                 const Vector<int> &moment_orders, SmallDataIO &file, 
                                 const int is_first_step);

    private:
        InflationConfig m_params;
        const int print_mode_functions = 0;
        std::string m_data_path;
        Real m_dt;
        Real m_cur_time;
        int m_restart_time;
        int m_first_step;
        const Real norm;

        std::string make_subdirectory(const std::string base, const std::string dir, 
                                      const int is_first_step) const;

        void assign_statistics_data(Vector<std::string> &header_storage, const std::string name, 
                                    Vector<Real> &data_storage, const Vector<Real> data, 
                                    const int component, const int num_comps,
                                    const Vector<int>::const_iterator itr, 
                                    const Vector<int>::const_iterator start, 
                                    const int is_first_step);

        void print_power_spectrum(const cMultiFab &field_array, SmallDataIO &power_spec_file, 
                                  const int component);
        Real calculate_field_moment_x(const MultiFab &field, const Vector<Real> mean, 
                                      const int moment, const int component);

        void extract_hs_and_R(MultiFab &hs, MultiFab &R, 
                              const MultiFab &state, const bool print_spec);
};

#include "InflationExtraction.impl.hpp"

#endif /* INFLATIONEXTRACTION_HPP_ */