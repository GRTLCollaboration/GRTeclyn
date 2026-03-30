/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef TENSORTESTS_HPP_
#define TENSORTESTS_HPP_

#include "InflationUtils.hpp"

using namespace amrex;

namespace TensorTests 
{
    // Test that the input tensor field (config space) is trace free (global)
    inline void Test_is_trace_free(MultiFab &field)
    {
        if (field.nComp() != 6)
        {
            Error("RandomField::Test_is_trace_free, input field is not a tensor field");
        }

        const auto &arrs = field.arrays();
        ParallelFor(field,
            [=] AMREX_GPU_DEVICE (int bx, int i, int j, int k)
            {
                IntVect iv{i, j, k};
                Real sum = 0.;

                for(int l=0; l<3; l++)
                {
                    sum += arrs[bx](i, j, k, InflationUtils::lut[l][l]);
                }

                if(std::abs(sum) > InflationUtils::tolerance)
                {
                    Print() << iv << ": " << sum << "\n";
                    Error("RandomField::Test_is_trace_free Trace-free test failed here.");
                }
            }
        );
    }

    // Test that the input vectors are orthonormal (local)
    inline void Test_vector_orthonorm(const IntVect iv, const Vector<Real> mhat, 
                                            const Vector<Real> nhat)
    {
        // Confirm basis vectors are orthonormal
        if (iv != IntVect{0, 0, 0})
        {
            Real dot1 = 0.;
            Real dot2 = 0.;
            Real cross1 = 0.;
            for(int l=0; l<3; l++)
            {
                dot1 += mhat[l] * mhat[l];
                dot2 += nhat[l] * nhat[l];
                cross1 += mhat[l] * nhat[l];
            }

            if(std::abs(dot1 - 1.) > InflationUtils::tolerance 
              || std::abs(dot2 - 1.) > InflationUtils::tolerance 
              || cross1 > InflationUtils::tolerance)
            {
                Print() << "Location: " << iv << "\n";
                Print() << "Dot products: " << dot1 << ", " << dot2 << "\n";
                Print() << "Cross product: " << cross1 << "\n";
                Print() << "Vectors: \n";
                for(int l=0; l<3; l++)
                {
                    Print() << l << ", " << mhat[l] << ", " << nhat[l] << "\n";
                }
                Error("RandomField::Test_vector_orthonorm: Basis vectors are not orthonormal here");
            }
        }
    }

    // Test that the input basis tensors, and their rotated counterparts, are orthonormal
    inline void Test_polarisation_tensor_orthonorm(const IntVect iv, const Tensor<2, Real> eplus,
                                                                    const Tensor<2, Real> ecross)
    {
        Vector<Real> conds(3, 0.);

        for (int l=0; l<3; l++) for (int p=0; p<3; p++)
        {
            conds[0] += eplus[l][p] * eplus[l][p];
            conds[1] += eplus[l][p] * ecross[l][p];
            conds[2] += ecross[l][p] * ecross[l][p];
        }

        if(iv != IntVect{0, 0, 0})
        {
            bool plc = (std::abs(conds[0] / 2. - 1.) > InflationUtils::tolerance 
                        || std::abs(conds[2] / 2. - 1.) > InflationUtils::tolerance);
            bool ppc = (std::abs(conds[1]) > InflationUtils::tolerance);
            if (plc || ppc)
            {
                Print() << "---------\nLocation: " << iv << "\n";
                Print() << "Dot products: " << conds[0] / 2. << ", " << conds[2] / 2. << "\n";
                Print() << "Cross product: " << conds[1] << "\n";
                Print() << "Base tensor components: \n";
                for (int l=0; l<3; l++) for (int p=0; p<3; p++)
                {
                    Print() << l << ", " << p << ": ";
                    Print() << eplus[l][p] << ", " << ecross[l][p] << "\n";
                }
                Error("RandomField::Test_polarisation_tensor_orthonorm: polarisation tensors are not orthonormal here");
            }
        }
    }
}

#endif /* TENSORTESTS_HPP_ */