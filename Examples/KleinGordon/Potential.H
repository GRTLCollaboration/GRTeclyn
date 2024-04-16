#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>

#include <cmath>

class Potential
{
  private:
  amrex::Real m_mass;
  
  public:
  Potential(const amrex::Real mass):m_mass(mass){};
  ~Potential(){};
  
  amrex::Real phi_sq(const amrex::Vector<amrex::Real> phi) const;
  amrex::Real sine_gordon(const amrex::Vector<amrex::Real> phi) const;

  
};
