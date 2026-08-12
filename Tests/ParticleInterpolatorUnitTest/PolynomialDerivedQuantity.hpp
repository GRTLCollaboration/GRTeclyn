/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef POLYNOMIALDERIVEDQUANTITY_HPP_
#define POLYNOMIALDERIVEDQUANTITY_HPP_

// AMReX includes
#include <AMReX_AmrLevel.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include "AMReX_Array.H"

// My attempt to make a derived variable (a polynomial). I followed what has
// been done for Constraints before.

class PolynomialDerivedQuantity
{
  public:

    static inline amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_center{
        AMREX_D_DECL(0., 0., 0.)}; // some random default here

    // read in the center which will be used to populating the polynomial
    // properly
    static void
    set_center(const std::array<amrex::Real, AMREX_SPACEDIM> &my_center)
    {
        for (int index = 0; index < AMREX_SPACEDIM; ++index)
        {
            m_center[index] = my_center[index];
        }
    }

    static inline const std::string name = "polynomial";

    static inline const amrex::Vector<std::string> var_names = {"pol"};

    // register with AMReX derive list
    AMREX_FORCE_INLINE static void set_up(int a_state_index)
    {
        const auto &desc_lst = amrex::AmrLevel::get_desc_lst();
        auto &derive_lst     = amrex::AmrLevel::get_derive_lst();

        derive_lst.add(name, amrex::IndexType::TheCellType(), 1, var_names,
                       compute_mf, amrex::DeriveRec::TheSameBox,
                       &amrex::cell_quartic_interp);

        // arguments: string (name), DescriptorList, int for state index,
        // starting comp, number of comps: we have only 1 polynomial so s_comp
        // (starting component) = 0.
        derive_lst.addComponent(name, desc_lst, a_state_index, 0, 1);
    }

    // Compute mf function, this is our DeriveFunc (bcrec is Boundary Condition
    // Records). As you can see some things are unused in this function, but we
    // need to provide them in the arguments as otherwise it does not compile :/
    static void compute_mf(amrex::MultiFab &out_mf, int dcomp, int ncomp,
                           const amrex::MultiFab &src_mf,
                           const amrex::Geometry &geom, amrex::Real /*time*/,
                           const int * /*bcrec*/, int /*level*/)
    {
        AMREX_ALWAYS_ASSERT(ncomp == 1);

        const auto plo = geom.ProbLoArray();
        const auto dx  = geom.CellSizeArray();

        const auto &out_arrays = out_mf.arrays();       // this is writable
        const auto &src_arrays = src_mf.const_arrays(); // this is read-only
        int ipoly              = dcomp;

        PolynomialDerivedQuantity polynomial(ipoly);

        auto center = m_center;

        amrex::ParallelFor(out_mf, out_mf.nGrowVect(),
                           [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                           {
                               // compute for given cell (i,j,k) in box number
                               // box_no
                               polynomial.compute(i, j, k, out_arrays[box_no],
                                                  src_arrays[box_no], plo, dx,
                                                  center);
                           });
    }

    // Compute function
    AMREX_GPU_DEVICE
    void
    compute(int i, int j, int k, const amrex::Array4<amrex::Real> &a_array,
            const amrex::Array4<amrex::Real const> &state,
            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> &plo,
            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> &dx,
            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> &center) const
    {
        // set up the coords
        amrex::Real x = plo[0] + (i + 0.5) * dx[0] - center[0];
        amrex::Real y = plo[1] + (j + 0.5) * dx[1] - center[1];
#if AMREX_SPACEDIM == 3
        amrex::Real z = plo[2] + (k + 0.5) * dx[2] - center[2];
#else
        amrex::Real z = 0.0;
#endif

        // write via cell data
        auto cell      = a_array.cellData(i, j, k);
        cell[m_c_poly] = 42.0 + x * x + y * y * z * z;
    }

    // Constructor
    PolynomialDerivedQuantity(int a_c_poly) : m_c_poly(a_c_poly) {}

  private:
    int m_c_poly; // destination comp
};

#endif /* POLYNOMIALDERIVEDQUANTITY_HPP_ */
