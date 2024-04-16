#ifndef INITIAL_CONDITIONS_H_
#define INITIAL_CONDITIONS_H_



#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <cmath>

class InitialConditions
{
 private:
  const amrex::Real m_alpha; 
  const amrex::Real m_beta; 

  const amrex::Real m_k_r; 
  const amrex::Real m_omega; 

  public:
  InitialConditions(const amrex::Real alpha, const amrex::Real k_r):m_alpha(alpha), m_beta(std::sqrt(1.0-alpha*alpha)), m_k_r(k_r), m_omega(k_r){};

  amrex::Real breather_solution(const amrex::Real x, const amrex::Real t) const;
  amrex::Real breather_solution(const amrex::Real x, const amrex::Real y, const amrex::Real z, const amrex::Real t) const;
  amrex::Real breather_solution_deriv(const amrex::Real x, const amrex::Real t) const;


  amrex::Real travelling_wave(const amrex::Real x, const amrex::Real y, const amrex::Real z, const amrex::Real t) const;
  amrex::Real travelling_wave_deriv(const amrex::Real x, const amrex::Real y, const amrex::Real z, const amrex::Real t) const;
 
};

#endif /* INITIAL_CONDITONS_H_ */
