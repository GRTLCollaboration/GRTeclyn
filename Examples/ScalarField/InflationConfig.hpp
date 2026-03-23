/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef INFLATIONCONFIG_HPP_
#define INFLATIONCONFIG_HPP_

#include "Cell.hpp"
#include "InitialScalarData.hpp"
#include "VarsTools.hpp"
#include "FilesystemTools.hpp"
#include "Potential.hpp"
#include "Tensor.hpp"
#include <fstream>
#include <random>

#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_FFT.H>
#include <AMReX_Random.H>
#include <AMReX_Print.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>

using namespace amrex;

struct InflationConfig
{
    // Basic initialisation flags
    int read_from_stoiic = 0;   //!< Whether to read spectrum from stoiic dparams.txt input
    int tensor_init = 0;        //!< Determines whether tensor perturbations are calculated
    int scalar_init = 0;  //!< Read in perturbations from STOIIC dparams
    int use_rand = 1;           //!< Choose whether to use random initial conditions
    int plot_int = 0;

    // Grid parameters
    Real L = 0;                   //!< Length of the box
    Real A = 0;                   //!< Amplitude factor (for basic tests)
    Real Mp = 1.;             //!< Energy scale of the problem
    Real alpha = 0.;          //!< Internal rotation angle in the +/x decomposition basis
    int N = 0;               //!< used to read in the private N variable
    int N_fine = 0;                 //!< Fine resolution to downsample from, 
                                //!< used for convergence testing

    // Initial condition options
    int random_seed = 3539263;  //!< Seed for random number generator
    int use_window = 0;         //!< Choose whether to use window function
    Real kstar = 0;               //!< window's cut-off mode, measured in units of 2pi/L
    Real Delta = 0;               //!< window's width, measured like L/Delta

    // Extraction parameters
    int calc_binned_power_spectrum = 0;   //!< Choose whether to extract the binned power spectrum
    int bin_number = 0;                       //!< How many bins to use (capped at N/2)
    int calc_higher_order_statistics = 0; //!< Choose whether to print higher-order statistics on the fields
    int num_orders = 0;                       //!< Number of moments to print (required by vector read-in)
    Vector<int> orders;                   //!< Moment orders to print for extracted fields
    int print_mode_functions = 0;

    // STOIIC read-in structures
    Vector<Real> init_k;                  //!< ks printed by STOIIC, at which Fourier-space fields are provided
    Vector<Vector<Real>> scalar_ps;       //!< Structure: four fields * two components, power spec values
    Vector<Vector<Real>> tensor_ps;       //!< Structure: two fields * two components, power spec values

    // Nyquist condition
    inline int flip_index(const int indx) 
    {
        AMREX_ASSERT(N > 0); 
        return std::abs(N - indx); 
    }

    // Nyquist condition and calculation of kmag
    inline int invert_index(const int indx) 
    { 
        AMREX_ASSERT(N > 0);
        return (int)(N/2 - std::abs(N/2 - indx)); 
    }

    // For calculation of polarisation tensors
    inline int invert_index_with_sign(const int indx) 
    { 
        AMREX_ASSERT(N > 0);
        if(indx <= N/2) { return indx; }
        else { return std::abs(N/2 - indx) - N/2; }
    }

    // Find the magnitude of the Fourier wavevector at this point
    inline Real get_kmag(IntVect iv)
    {
        AMREX_ASSERT(L > 0);
        const int i = iv[0];
        const int j = invert_index(iv[1]);
        const int k = invert_index(iv[2]);
        return std::sqrt(i*i + j*j + k*k) * 2. * M_PI / L;
    }

    inline Vector<Real> InflationConfig::calculate_basis_vector(const IntVect iv, 
                                                        const int which_vector);

    // Applies above Nyquist conditions to a given MF
    inline void InflationConfig::apply_nyquist_conditions(cMultiFab &field);
}