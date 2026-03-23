/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef RANDOMFIELDINIT_HPP_
#define RANDOMFIELDINIT_HPP_

#include "InflationConfig.hpp"
#include "TensorTests.hpp"
using namespace amrex;

class RandomFieldInit
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

        const Real& H0;
        static Real calc_H0(Real G, Real Pi, Real V)
        {
            return sqrt((8. * M_PI * G / 3.)*(0.5*pow(Pi, 2.) + V));
        }

    public:
                // Constructor used when initialising stochastic fields
        RandomFieldInit(const InflationConfig a_config, 
                    const InitialBackgroundData::params_t bkgd_params, 
                    const Potential::params_t potential_params)
                : m_params(a_config), tol(1.e-12), 
                  norm(pow(1./a_params.L, 3.))
        {
            // Compute background potential
            double V, dV;
            Potential potential(potential_params);
            switch (potential_params.type)
            {
                case 1:
                    potential.quadratic(V, dV, bkgd_params.phi0);
                case 9:
                    potential.monodromy(V, dV, bkgd_params.phi0);
                case 8:
                    potential.USR(V, dV, bkgd_params.phi0);
                case 10:
                    potential.punctuated(V, dV, bkgd_params.phi0);
                default:
                    Error("RandomFieldInit::RandomFieldInit, potential type not provided");
            }

            // Compute initial Hubble parameter
            H0 = calc_H0(bkgd_params.G_Newton, bkgd_params.Pi, V);
        }

        void init(amrex::MultiFab &state);

    private:
        void make_random_draws(MultiFab &rand_fab, Box &domain, const int seed);
        GpuComplex<Real> calculate_mode_function(const double km, const int spec_indx);
        GpuComplex<Real> find_in_stoiic(const double km, const int field_indx, 
                                        std::string field_type);
        GpuComplex<Real> calculate_random_field(const IntVect iv, const int field_index, 
                                                const Real rand_amp, const Real rand_phase, 
                                                std::string field_type);
};

#include "RandomFieldInit.impl.hpp"

#endif /* RANDOMFIELDINIT_HPP_ */