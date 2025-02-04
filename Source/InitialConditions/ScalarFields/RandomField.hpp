/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef RANDOMFIELD_HPP_
#define RANDOMFIELD_HPP_

#include "Cell.hpp"
#include "InitialScalarData.hpp"
#include "VarsTools.hpp"
#include "FilesystemTools.hpp"
#include <fstream>

#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_FFT.H>
#include <AMReX_Random.H>
#include <AMReX_Print.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>

using namespace amrex;

//! Class to create a Gaussian random field, 
//! originally created for 2 massless tensor polarisation fields
//! but will be extended to N IID fields with given masses.
class RandomField 
{
    public:
        //! A structure for storing parameters essential to this class
        struct params_t 
        {
            int num_scalar_fields;      //!< Number of fields to generate
            int calc_tensor_field;      //!< Determines whether tensor perturbations are calculated
            int use_rand = 1;           //!< Choose whether to use random initial conditions
            int random_seed = 3539263;  //!< Seed for random number generator

            double L;                   //!< Length of the box
            double A;                   //!< Amplitude factor (for basic tests)
            double Mp = 1.;             //!< Energy scale of the problem

            int N_readin;               //!< used to read in the private N variable
            int N_fine;                 //!< Fine resolution to downsample from, 
                                        //!< used for convergence testing

            int use_window = 0;         //!< Choose whether to use window function
            double kstar;               //!< window's cut-off mode, measured in units of 2pi/L
            double Delta;               //!< window's width, measured like L/Delta

            int calc_binned_power_spectrum = 0;   //!< Choose whether to extract the binned power spectrum
            int bin_number = N_readin/2;          //!< How many bins to use (capped at N/2)
            int calc_config_space_mode_fns = 0;   //!< Choose whether to print the fields in configuration space
            int calc_higher_order_statistics = 0; //!< Choose whether to print higher-order statistics on the fields
            int num_orders;
            Vector<int> orders;                   //!< Moment orders to print for extracted fields
        };

        RandomField(params_t a_params, InitialBackgroundData::params_t a_background_params)
                : m_params(a_params), m_background_params(a_background_params)
        {
            // Set protected class parameters
            N = m_params.N_readin;
            norm = m_params.A * pow(2. * M_PI/m_params.L, 3.); // Physical FFT normalisation

            // Look-up table 
            // Used to construct polarisation basis tensors
            lut[0][0] = 0;
            lut[0][1] = 1;
            lut[0][2] = 2;
            lut[1][0] = 1;
            lut[1][1] = 3;
            lut[1][2] = 4;
            lut[2][0] = 2;
            lut[2][1] = 4;
            lut[2][2] = 5;
        }

        void init(amrex::MultiFab &state);
        void extract(MultiFab &state, std::string data_path, Real dt, 
                        Real cur_time, int restart_time, int first_step);
        
    private:
        int N;
        int lut[3][3];
        double norm;

        int flip_index(int indx);
        int invert_index(int indx);
        int invert_index_with_sign(int indx);
        std::string make_subdirectory(std::string base, std::string dir, int is_first_step);
        void assign_statistics_data(Vector<std::string> &header_storage, const std::string name, 
                            Vector<Real> &data_storage, const Array1D<Real, 0, 1> data, 
                            const int component, const auto itr, const auto start, const int is_first_step);

        GpuComplex<Real> calculate_mode_function(double km, std::string spec_type);
        GpuComplex<Real> calculate_random_field(int I, int J, int k, std::string spectrum_type, 
                                                                Real rand_amp, Real rand_phase);
        Vector<Real> calculate_basis_vector(int i, int j, int k, int which_vector);
        GpuComplex<Real> calculate_tensor_initial_conditions(int I, int J, int k, int l, int p, 
                            GpuComplex<Real> plus_field, GpuComplex<Real> cross_field);
        bool is_independent_draw(IntVect iv);
        void apply_nyquist_conditions(cMultiFab &field, BaseFab<GpuComplex<Real>> &plane1, BaseFab<GpuComplex<Real>> &plane2);
        bool is_ghost_index(IntVect vector);

        Real find_field_moment_x(MultiFab &field, Array1D<Real, 0, 1> mean, 
                                                int moment, int component);
        void print_tensor_moment(MultiFab &field, const Vector<int> &moment_orders, 
                                    SmallDataIO &file, const int is_first_step);
        void print_power_spectrum(cMultiFab &field_array, SmallDataIO &power_spec_file, int component);

    protected:
        const params_t m_params;
        const InitialBackgroundData::params_t m_background_params;
        const std::string m_spec_type;
};

#include "RandomField.impl.hpp"

#endif /* RANDOMFIELD_HPP_ */
