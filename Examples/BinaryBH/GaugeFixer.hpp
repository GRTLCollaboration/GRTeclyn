/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef GAUGEFIXER_HPP_
#define GAUGEFIXER_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "StateVariables.hpp"
#include "Tensor.hpp"

// This class fixes the gauge far away
class GaugeFixer
{
  public:

    // Constructor
    GaugeFixer(double a_dx, std::array<double, AMREX_SPACEDIM> a_center) : 
	    m_dx(a_dx), m_center(a_center)
    {}

    amrex::Real m_dx;
    std::array<double, AMREX_SPACEDIM> m_center;

    // Compute function
    AMREX_GPU_DEVICE void
    operator()(int ix, int iy, int iz,
               const amrex::Array4<amrex::Real> &rhs_state) const
    {
        Coordinates coords(amrex::IntVect(ix, iy, iz), m_dx);
    
	const amrex::CellData<amrex::Real> &rhs_cell_data =
            rhs_state.cellData(ix, iy, iz);

	amrex::Real r = std::sqrt(std::pow(coords.x - m_center[0], 2) +
                              std::pow(coords.y - m_center[1], 2) +
                              std::pow(coords.z - m_center[2], 2));

        r = std::max(r, 1e-6);
        const double R = 500.0;
        const double f_of_r = (R * R / (r * r + R * R));

        FOR(i)
        {
            rhs_cell_data[c_B1 + i] *= f_of_r ;
        }
    }
};

#endif /* GAUGEFIXER_HPP_ */
