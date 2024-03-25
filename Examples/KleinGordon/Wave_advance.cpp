#include "KleinGordonLevel.hpp"

using namespace amrex;

Real
KleinGordonLevel::advance (Real time, Real dt, int iteration, int ncycle)
{
    // At the beginning of step, we make the new data from previous step the
    // old data of this step.
    for (int k = 0; k < NUM_STATE_TYPE; ++k) {
        state[k].allocOldData();
        state[k].swapTimeLevels(dt);
    }

    AmrLevel::RK(rk_order, State_Type, time, dt, iteration, ncycle,
                 [&] (int /*stage*/, MultiFab& dSdt, MultiFab const& S,
                      Real /*t*/, Real /*dtsub*/)
    {
        computeRHS(dSdt, S);
    });

    return dt;
}

// Real 
// AmrLevelWave::Potential(Real phi2)
// {
//   Real V = 0;
//   //  Real scalar_mass = 0;

//   V = 0.5*scalar_mass*scalar_mass*phi2;
//   //  V = std::sin(phi2);

//   return V;
// }


void
KleinGordonLevel::computeRHS (MultiFab& dSdt, MultiFab const& S)
{
    const auto dxinv = Geom().InvCellSizeArray();
    AMREX_D_TERM(Real dx2inv = dxinv[0]*dxinv[0];,
                 Real dy2inv = dxinv[1]*dxinv[1];,
                 Real dz2inv = dxinv[2]*dxinv[2]);
    auto const& sa = S.const_arrays();
    auto const& sdot = dSdt.arrays();

    Potential my_potential(scalar_mass);
	  
    amrex::ParallelFor(S,
    [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k) noexcept
    {
      auto const& s = sa[bi];
      auto const& f = sdot[bi];



      //      Real phi2 = std::pow(s(i,j,k,0),2)+std::pow(s(i,j,k,2),2);
      amrex::Vector<amrex::Real> phi;


      for (int n = 0; n < nfields; n++)
	{
	  phi.push_back(s(i,j,k,2*n));

	  f(i,j,k,2*n) = s(i,j,k,2*n+1);

	  AMREX_D_TERM(Real lapx = dx2inv*(-2.5*s(i,j,k,2*n) + (4./3.)*(s(i-1,j,k,2*n)+s(i+1,j,k,2*n))
					   -                (1./12.)*(s(i-2,j,k,2*n)+s(i+2,j,k,2*n)));,
		       Real lapy = dy2inv*(-2.5*s(i,j,k,2*n) + (4./3.)*(s(i,j-1,k,2*n)+s(i,j+1,k,2*n))
					   -                (1./12.)*(s(i,j-2,k,2*n)+s(i,j+2,k,2*n)));,
		       Real lapz = dz2inv*(-2.5*s(i,j,k,2*n) + (4./3.)*(s(i,j,k-1,2*n)+s(i,j,k+1,2*n))
					   -                (1./12.)*(s(i,j,k-2,2*n)+s(i,j,k+2,2*n))));

	  f(i,j,k,2*n+1) = AMREX_D_TERM(lapx, +lapy, +lapz);

	  f(i,j,k,2*n+1) -= std::sin(s(i,j,k,2*n));

	}


    });
    Gpu::streamSynchronize();
}
