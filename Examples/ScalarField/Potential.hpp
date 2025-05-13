/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef POTENTIAL_HPP_
#define POTENTIAL_HPP_

#include "simd.hpp"

class Potential
{
  public:
    struct params_t
    {
        double scalar_mass;

	double phi0;
	double df0;
	double Mp; 

	double width;	
	double location;
	double amplitude;
	double wavelength;
    };

  private:
    params_t m_params;

  public:
    //! The constructor
    Potential(params_t a_params) : m_params(a_params) 
    {
	m_params.scalar_mass *= m_params.Mp;
	m_params.phi0 *= m_params.Mp;
	m_params.df0 *= std::pow(m_params.Mp, 2.);

	m_params.width = 2. * m_params.df0 / (9.e-6) / m_params.Mp;
	m_params.location = m_params.phi0 + 0.5 * m_params.width;
	m_params.amplitude *= std::pow(m_params.Mp, 4.);
	m_params.wavelength *= m_params.Mp;
    }

    //! Set the potential function for the scalar field here
    template <class data_t, template <typename> class vars_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    compute_potential(data_t &V_of_phi, data_t &dVdphi,
                      const vars_t<data_t> &vars) const
    {
	// The potential value at phi
	// Monodromy model, loosely based off the one used in STOIIC_GR
	double argument = (vars.phi - m_params.location)/m_params.wavelength;
	double displaced_argument = (m_params.location - vars.phi + m_params.width)/m_params.wavelength;
	
	double envelope = 0.25 * (1. + tanh(argument)) * (1. + tanh(displaced_argument));
	double oscillation = cos(argument) - 1.; 	

	V_of_phi = 0.5 * pow(m_params.scalar_mass * vars.phi, 2.0);
	V_of_phi += m_params.amplitude * (envelope * oscillation);

        // The potential gradient at phi
	double d_envelope = 0.25/m_params.wavelength * 
				((1. + tanh(argument)) * (std::pow(tanh(displaced_argument), 2.) - 1.)
			       + (1. + tanh(displaced_argument)) * (1. - std::pow(tanh(argument), 2.)));
	double d_oscillation = -sin(argument)/m_params.wavelength;

        dVdphi = pow(m_params.scalar_mass, 2.0) * vars.phi;
	dVdphi += m_params.amplitude * (envelope * d_oscillation + d_envelope * oscillation);
    }
};

#endif /* POTENTIAL_HPP_ */
