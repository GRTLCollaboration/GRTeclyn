/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef INFLATIONEXTRACTION_HPP_
#define INFLATIONEXTRACTION_HPP_

#include "InflationConfig.hpp"
using namespace amrex;

class InflationExtraction : RandomFieldInit
{
    protected:
        const Real tol;
        const Real norm;

        // Look-up table 
        // Used to construct polarisation basis tensors
        int temp_lut[3][3];
        temp_lut[0][0] = 0;
        temp_lut[0][1] = 1;
        temp_lut[0][2] = 2;
        temp_lut[1][0] = 1;
        temp_lut[1][1] = 3;
        temp_lut[1][2] = 4;
        temp_ut[2][0] = 2;
        temp_lut[2][1] = 4;
        temp_lut[2][2] = 5;
        const int lut[3][3] = temp_lut;

    public:
        // Names of diagnostic variables
        static inline const Vector<std::string> var_names = {"hplus", "hcross"};

        // Constructor used in extraction of diagnostics
        RandomField(InflationConfig::InflationConfig a_params)
                : m_params(a_params), tol(1.e-12), 
                  norm(pow(1./a_params.L, 3.))       
        {}

        void derive(const MultiFab &source, MultiFab &out, int dcomp);
        void extract(const MultiFab &state, const std::string data_path, const Real dt,  
                     const Real cur_time, const int restart_time, const int first_step);
        Vector<Real> print_moment(MultiFab &field, const Vector<std::string> names,  
                                 const Vector<int> &moment_orders, SmallDataIO &file, 
                                 const int is_first_step);

    private:
        std::string make_subdirectory(const std::string base, const std::string dir, 
            const int is_first_step);
        void assign_statistics_data(Vector<std::string> &header_storage, const std::string name, 
                                    Vector<Real> &data_storage, const Vector<Real> data, 
                                    const int component, const int num_comps,
                                    const Vector<int>::const_iterator itr, 
                                    const Vector<int>::const_iterator start, 
                                    const int is_first_step);

        void print_power_spectrum(cMultiFab &field_array, SmallDataIO &power_spec_file, 
                                    const int component);
        Real find_field_moment_x(MultiFab &field, const Vector<Real> mean, 
                                 const int moment, const int component);
}