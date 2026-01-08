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
		int type;
		double param1, param2, param3, param4, param5;

		// Monodromy parameters
		double scalar_mass = 0.;
		double location;
		double width;	
		double amplitude;
		double period;

		// USR (Prokopec) parameters
		double Lambda;
		double v = 0.;
    };

  private:
    params_t m_params;

  public:
    //! The constructor
    Potential(params_t a_params) : m_params(a_params) 
	{
		switch (m_params.type)
		{
			case 1:
				m_params.scalar_mass = m_params.param1;

			case 8:
				m_params.scalar_mass = m_params.param1;
				m_params.location = m_params.param2;
				m_params.width = m_params.param3;
				m_params.amplitude = m_params.param4;
				m_params.period = m_params.param5;

			case 9:
				m_params.Lambda = m_params.param1;
				m_params.v = m_params.param2;
		}
	}

	// Classic quadratic potenital
	template <class data_t>
	AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
	quadratic(data_t &V, data_t &dV, const data_t &phi) const
	{
		if (m_params.scalar_mass == 0)
		{
			amrex::Error("Potential::quadratic, Scalar mass is un-initialised.");
		}

		V = 0.5 * std::pow(m_params.scalar_mass * phi, 2.);
		dV = std::pow(m_params.scalar_mass, 2.) * phi;
	}

	// Monodromy potential, as used in STOIIC and also in arXiv:2403.12811
	template <class data_t>
	AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
	monodromy(data_t &V, data_t &dV, const data_t &phi) const
	{
		if (m_params.scalar_mass == 0)
		{
			amrex::Error("Potential::monodromy, Scalar mass is un-initialised.");
		}

		// Calculate V
		double argument = (phi - m_params.location)/m_params.period;
		double displaced_argument = (m_params.location - phi + m_params.width)/m_params.period;
		
		double envelope = 0.25 * (1. + tanh(argument)) * (1. + tanh(displaced_argument));
		double oscillation = cos(argument) - 1.; 	

		V = 0.5 * pow(m_params.scalar_mass * phi, 2.0);
		V += m_params.amplitude * (envelope * oscillation);

		// Calculate dV
		double d_envelope = 0.25/m_params.period * 
					((1. + tanh(argument)) * (std::pow(tanh(displaced_argument), 2.) - 1.)
					+ (1. + tanh(displaced_argument)) * (1. - std::pow(tanh(argument), 2.)));
		double d_oscillation = -sin(argument)/m_params.period;

		dV = pow(m_params.scalar_mass, 2.0) * phi;
		dV += m_params.amplitude * (envelope * d_oscillation + d_envelope * oscillation);
	}

	// Prokopec USR model, from arXiv:2507.04114
	template <class data_t>
	AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
	USR(data_t &V, data_t &dV, const data_t &phi) const
	{
		if (m_params.v == 0)
		{
			amrex::Error("Potential::USR, USR parameter v is un-initialised.");
		}

		// Calculate V
		double fraction = (3. * pow(phi, 2.) 
						 + 2. * sqrt(2.) * phi * m_params.v 
						 + 6. * pow(m_params.v, 2.));

		fraction /= pow(3. * pow(phi, 2.) + 2. * pow(m_params.v, 2.), 2.);
		V = m_params.Lambda * pow(m_params.v, 4.) * pow(phi, 2.) * fraction / 3.;

		// Calculate dV
		fraction = ((2. * m_params.v + sqrt(2.) * phi) 
				  * (pow(phi, 2.) - 2. * pow(m_params.v, 2.)));
		fraction /= pow(2. * pow(m_params.v, 2.) + 3. * pow(phi, 2.), 3.);
		dV = - 2. * m_params.Lambda * pow(m_params.v, 5.) * phi * fraction;
	}

    //! Set the potential function for the scalar field here
    template <class data_t, template <typename> class vars_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    compute_potential(data_t &V_of_phi, data_t &dVdphi,
                      const vars_t<data_t> &vars) const
    {
		if(m_params.type == 1)
		{
			quadratic<data_t>(V_of_phi, dVdphi, vars.phi);
		}
		else if (m_params.type == 8)
		{
			USR(V_of_phi, dVdphi, vars.phi);
		}
		else if (m_params.type == 9)
		{
			monodromy(V_of_phi, dVdphi, vars.phi);
		}
		else
		{
			amrex::Print() << m_params.type << "\n";
			amrex::Error("Potential::compute_potential, requested potential type is not supported.");
		}
    
		/*amrex::Print().SetPrecision(15) << "V: " << V_of_phi << "\n";
		amrex::Print().SetPrecision(15) << "dV: " << dVdphi << "\n";
		amrex::Error();*/
    }
};

#endif /* POTENTIAL_HPP_ */
