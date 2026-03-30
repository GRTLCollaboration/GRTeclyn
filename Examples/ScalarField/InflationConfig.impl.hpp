/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */


#if !defined(INFLATIONCONFIG_HPP_)
#error "This file should only be included via InflationConfig.hpp"
#endif

#ifndef INFLATIONCONFIG_IMPL_HPP_
#define INFLATIONCONFIG_IMPL_HPP_

// Calculates basis vectors required for polarisation tensors
inline Vector<Real> InflationConfig::calculate_basis_vector(const IntVect iv, const int which_vector)
{
    AMREX_ASSERT(norm > 0);

    // Hermitian symmetry inversion on j and k, with sign on the last two indices.
    // (!!) The FT implemented in AMReX symmetrises across the i index.
    const Real i = static_cast<Real>(iv[0]);
    const Real j = static_cast<Real>(invert_index_with_sign(iv[1]));
    const Real k = static_cast<Real>(invert_index_with_sign(iv[2]));

    Vector<Real> mhat(3, 0.);
    Vector<Real> nhat(3, 0.);

    // Skip the 0 mode, as tensors have no average
    if (iv == IntVect{0, 0, 0}) { ; }

    else if (i > 0.) 
    {
        if (k == 0. && j == 0.) 
        { 
            mhat = Vector<Real>{0., 1., 0.};
            nhat = Vector<Real>{0., 0., 1.}; 
        }

        else 
        { 
            Real norm = sqrt((i*i + j*j) * (i*i + j*j + k*k));
            mhat = Vector<Real>{j/sqrt(i*i + j*j), -i/sqrt(i*i + j*j), 0.}; 
            nhat = Vector<Real>{(i*k) / norm, (j*k) / norm, -(i*i + j*j) / norm}; 
        }
    }

    else if (std::abs(j) > 0) 
    { 
        if(k == 0.)
        {
            mhat = Vector<Real>{0., 0., 1.};
            nhat = Vector<Real>{1., 0., 0.};
        }

        else
        {
            mhat = Vector<Real>{-1., 0., 0.};
            nhat = Vector<Real>{0., -k / sqrt(j*j + k*k), j / sqrt(j*j + k*k)};
        }
    }

    else if (std::abs(k) > 0) 
    { 
        mhat = Vector<Real>{1., 0., 0.};
        nhat = Vector<Real>{0., 1., 0.};
    }

    else 
    {
        Error("RandomField::calculate_polarisation_tensors Part of Fourier grid not covered.");
    }

    if (alpha != 0)
    {
        Real a = alpha * (M_PI) / 180.;
        Vector<Real> mp(3), np(3);
        for(int l=0; l<3; l++)
        {
            mp[l] = cos(a) * mhat[l] + sin(a) * nhat[l];
            np[l] = -sin(a) * mhat[l] + cos(a) * nhat[l];
        }

        mhat = mp;
        nhat = np;
    }

    if(which_vector == 0) { return mhat; }
    else if(which_vector == 1) { return nhat; }
    else 
    { 
        Error("RandomField::calculate_basis_vector Incompatable vector type."); 
        return Vector<Real>{0,0,0}; 
    }
}

// Applies above Nyquist conditions to a given MF
inline void InflationConfig::apply_nyquist_conditions(cMultiFab &field)
{
    AMREX_ASSERT(N > 0);
    
    int nc = field.nComp();
    for (MFIter mfi(field); mfi.isValid(); ++mfi) 
    {
        // The geometry for this MPI rank
        const Box& bx = mfi.fabbox();
        Array4<GpuComplex<Real>> const& field_ptr = field.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            IntVect iv = {i, j, k};

            if ((i == 0 || i == N/2) && (j == 0 || j == N/2) && (k == 0 || k == N/2))
            {
                for(int comp = 0; comp < nc; comp++)
                {
                    GpuComplex<Real> temp(field_ptr(i, j, k, comp).real(), 0.);
                    field_ptr(i, j, k, comp) = temp;
                }
            }

            else if (i==0 || i==N/2) 
            {
                if((k > N/2 && j == N/2) || (k == 0 && j > N/2) ||
                    (k > N/2 && j == 0) || (k == N/2 && j > N/2))
                {
                    for(int comp = 0; comp < nc; comp++) 
                    {
                        GpuComplex<Real> temp(field_ptr(i, invert_index(j), invert_index(k), comp).real(), 
                                                -field_ptr(i, invert_index(j), invert_index(k), comp).imag());
                        field_ptr(i, j, k, comp) = temp;
                    }
                }
                
                else if(j > N/2)
                {
                    for(int comp = 0; comp < nc; comp++) 
                    {
                        GpuComplex<Real> temp(field_ptr(i, invert_index(j), flip_index(k), comp).real(), 
                                                -field_ptr(i, invert_index(j), flip_index(k), comp).imag());
                        field_ptr(i, j, k, comp) = temp;
                    }
                }
            }
        });
    }
}

#endif /* INFLATIONCONFIG_IMPL_HPP */