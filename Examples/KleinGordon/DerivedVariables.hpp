#ifndef DERIVEDVARIABLES_HPP
#define DERIVEDVARIABLES_HPP

#include <AMReX_BLFort.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>

#include "VarsTools.hpp"

#include "InitialConditions.hpp"

void calc_analytic_solution(amrex::MultiFab &mf_out, int dcomp, int /*numcomp*/,
                            const amrex::MultiFab & /*mf_in*/,
                            const amrex::Geometry &geom, const amrex::Real time,
                            const int * /*bcomp*/, int /*scomp*/);

void calc_sine_gordon_1d_analytic_solution(amrex::MultiFab &mf_out, int dcomp,
                                           const amrex::Geometry &geom,
                                           const amrex::Real time);

void calc_sine_gordon_3d_analytic_solution(amrex::MultiFab &mf_out, int dcomp,
                                           const amrex::Geometry &geom,
                                           const amrex::Real time);

void calc_wave_analytic_solution(amrex::MultiFab &mf_out, int dcomp,
                                 const amrex::Geometry &geom,
                                 const amrex::Real time);

#endif // DERIVEDVARIABLES_HPP
