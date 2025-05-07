/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef INITIALBACKGROUNDDATA_HPP_
#define INITIALBACKGROUNDDATA_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "MatterCCZ4RHS.hpp"
#include "ScalarField.hpp"
#include "StateVariables.hpp" //This files needs NUM_VARS - total no. components
#include "Tensor.hpp"
#include "VarsTools.hpp"
#include "simd.hpp"

//#include "MayDay.H"
//#include <fstream>

class InitialBackgroundData
{
	public:
		struct params_t
		{
			double phi0; //!< Amplitude of k=0 mode of initial SF
			double Pi0;  //!< Amplitude of initial SF velocity
			double m;    //!< SF mass
			double Mp = 1.;    //!< Energy scale [Mp]
		};

		InitialBackgroundData(params_t a_params)
			: m_params(a_params)
		{
		}

		template <class data_t> 
		AMREX_GPU_DEVICE AMREX_FORCE_INLINE void 
		compute(int i, int j, int k, const amrex::Array4<data_t> &cell) const
		{
			MatterCCZ4RHS<ScalarField<>>::Vars<data_t> vars;
        		VarsTools::assign(vars, 0.); // Set only the non-zero components below

        		// start with unit lapse and flat metric (must be relaxed for chi)
        		vars.lapse = 1.0;
        		vars.chi   = 1.0;

			FOR(index)
				vars.shift[index] = 0.;
        		// conformal metric is flat
        		FOR (index)
            			vars.h[index][index] = 1.;

			const double phi = m_params.phi0 * m_params.Mp;
			const double Pi = m_params.Pi0 * std::pow(m_params.Mp, 2.);
			const double H0 = sqrt((8. * M_PI/3./pow(m_params.Mp, 2.))*0.5*(pow(Pi, 2.) 
						+ pow(m_params.m * m_params.Mp * phi, 2.0)));

			vars.phi = phi;
			vars.Pi = Pi;
			vars.K = -3.*H0;

			store_vars(cell.cellData(i, j, k), vars);
		}
	protected:
		const params_t m_params;

};


#endif /* INITIALBACKGROUNDDATA_HPP_ */
