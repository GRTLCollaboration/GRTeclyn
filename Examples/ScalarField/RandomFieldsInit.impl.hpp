/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */


#if !defined(INFLATIONEXTRACTION_HPP_)
#error "This file should only be included via RandomField.hpp"
#endif

#ifndef INFLATIONEXTRACTION_IMPL_HPP_
#define INFLATIONEXTRACTION_IMPL_HPP_


// Main initialisation routine
inline void RandomFieldInit::init(amrex::MultiFab &state)
{
    BL_PROFILE("RandomFieldInit::init");

    // Derive MultiFab ingredients from state (configuration space)
    BoxArray sba = state.boxArray();
    DistributionMapping sdm = state.DistributionMap();

    // Set up the problem domain in Fourier space
    // And impose that MPI ranks only slice along the i index (for Nyquist conditions)
    IntVect domain_low(0, 0, 0);
    IntVect k_domain_high(N/2, N-1, N-1);
    Box k_domain(domain_low, k_domain_high);
    Array< bool, AMREX_SPACEDIM > const &slicing{true, false, false};
    BoxArray kba = decompose(k_domain, ParallelContext::NProcsAll(), slicing);
    DistributionMapping kdm(kba);

    // Set up the MFs to store the in/out data sets
    cMultiFab hs_k(kba, kdm, 2, 0);
    cMultiFab As_k(kba, kdm, 2, 0);
    cMultiFab hij_k(kba, kdm, 6, 0);
    cMultiFab Aij_k(kba, kdm, 6, 0);

    MultiFab hij_x(sba, sdm, 6, 0);
    MultiFab Aij_x(sba, sdm, 6, 0);

    cMultiFab scalar_fields_k(kba, kdm, 4, 0);
    MultiFab scalar_fields_x(sba, sdm, 4, 0);

    hs_k.setVal(0.0);
    As_k.setVal(0.0);
    hij_k.setVal(0.0);
    Aij_k.setVal(0.0);
    hij_x.setVal(0.0);
    Aij_x.setVal(0.0);
    scalar_fields_k.setVal(0.0);
    scalar_fields_x.setVal(0.0);

    // Construct the Fourier transform
    IntVect x_domain_high(N-1, N-1, N-1);
    Box x_domain(domain_low, x_domain_high);
    FFT::R2C<Real> tensor_fft(x_domain, FFT::Info().setBatchSize(hij_k.nComp()));
    FFT::R2C<Real> scalar_fft(x_domain, FFT::Info().setBatchSize(scalar_fields_k.nComp()));

    // Construct MFs to hold the random number draws
    MultiFab random_draws(kba, kdm, 6, 0);
    make_random_draws(random_draws, k_domain, m_params.random_seed);
    MultiFab tensor_draws(random_draws, amrex::make_alias, 0, 4);
    MultiFab scalar_draws(random_draws, amrex::make_alias, 4, 2);

    const auto &tensor_draw_arrs = tensor_draws.arrays();
    const auto &scalar_draw_arrs = tensor_draws.arrays();
    const auto &hs_arrs = hs_k.arrays();
    const auto &hij_x_arrs = hij_k.arrays();
    const auto &As_arrs = As_k.arrays();
    const auto &Aij_x_arrs = Aij_k.arrays();
    const auto &scalar_field_arrs = scalar_fields_k.arrays();

    Print() << "Starting initial condition generation/read in...\n";
    
    ParallelFor(hs_k,
        [=] AMREX_GPU_DEVICE (int bx, int i, int j, int k)
    {
        IntVect iv = {i, j, k};
        if(m_params.scalar_init)
        {
            for(int f=0; f<4; f++)
            {
                Real draw1 = scalar_draw_arrs[bx](i, j, k, 0);
                Real draw2 = scalar_draw_arrs[bx](i, j, k, 1);

                scalar_fields_arrs[bx](i, j, k, f) = calculate_random_field(iv, f, draw1, draw2, "scalar");
            }
        }

        if(m_params.tensor_init)
        {
            // Find the mode function realisation
            for(int p=0; p<2; p++)
            {
                Real draw1 = tensor_draw_arrs[bx](i, j, k, 2*p);
                Real draw2 = tensor_draw_arrs[bx](i, j, k, 2*p+1);

                hs_arrs[bx](i, j, k, p) = calculate_random_field(iv, 0, draw1, draw2, "tensor");
                As_arrs[bx](i, j, k, p) = calculate_random_field(iv, 1, draw1, draw2, "tensor");
            }
            
            // Find basis vectors
            Vector<Real> mhat = m_params.calculate_basis_vector(iv, 0);;
            Vector<Real> nhat = m_params.calculate_basis_vector(iv, 1);
            TensorTests::Test_vector_orthonorm(iv, mhat, nhat);

            // Construct polarisation tensors from basis vectors
            Tensor<2, Real> eplus, ecross; 

            // Find basis tensors and initial tensor realisation
            for (int l=0; l<3; l++) for (int p=0; p<3; p++)
            {
                // Assemble the polarisation tensors
                eplus[l][p] = mhat[l]*mhat[p] - nhat[l]*nhat[p];
                ecross[l][p] = mhat[l]*nhat[p] + nhat[l]*mhat[p];

                hij_x_arrs[bx](i, j, k, lut[l][p]) = hs_arrs[bx](i, j, k, 0) * eplus[l][p] + hs_arrs[bx](i, j, k, 1) * ecross[l][p];
                Aij_x_arrs[bx](i, j, k, lut[l][p]) = As_arrs[bx](i, j, k, 0) * eplus[l][p] + As_arrs[bx](i, j, k, 1) * ecross[l][p];
            }

            TensorTests::Test_polarisation_tensor_orthonorm(iv, eplus, ecross);
        }
    });

    // Apply the DC and Nyquist symmetry conditions
    m_params.apply_nyquist_conditions(hs_k);
    m_params.apply_nyquist_conditions(hij_k);
    m_params.apply_nyquist_conditions(Aij_k);
    m_params.apply_nyquist_conditions(scalar_fields_k);

    // Do the Fourier transform
    tensor_fft.backward(hij_k, hij_x);
    tensor_fft.backward(Aij_k, Aij_x);
    scalar_fft.backward(scalar_fields_k, scalar_fields_x);

    // Apply normalisation into physical units
    hij_x.mult(norm);
    Aij_x.mult(norm);
    scalar_fields_x.mult(norm);

    // Test that the resuling tensor perturbation field is trace-free
    TensorTests::Test_is_trace_free(hij_x);

    // Convert to BSSN variables using the BSSN-CPT dictionary
    for (int l=0; l<3; l++) { hij_x.plus(1., lut[l][l], 1); }
    Aij_x.mult(-0.5);

    const auto &state_arrs = state.arrays();
    const auto &hij_x_arrs = hij_x.arrays();
    const auto &Aij_x_arrs = Aij_x.arrays();
    const autp &scalar_field_x_arrs = scalar_fields_x.arrays();

    ParallelFor(hij_x,
        [=] AMREX_GPU_DEVICE (int bx, int i, int j, int k) noexcept
    {
        const IntVect iv{i, j, k};
        // Add scalar perturbations to the existing background values
        if(m_params.scalar_init)
        {
            state_arrs[bx](iv, c_phi) += scalar_field_x_arrs[bx](i, j, k, 0);
            state_arrs[bx](iv, c_Pi) += scalar_field_x_arrs[bx](i, j, k, 1);
            state_arrs[bx](iv, c_chi) += scalar_field_x_arrs[bx](i, j, k, 2);
            state_arrs[bx](iv, c_K) += scalar_field_x_arrs[bx](i, j, k, 3);
        }

        // Set the entire tensor object here
        if(m_params.tensor_init)
        {
            state_arrs[bx](iv, c_h11) = hij_x_arrs[bx](i, j, k, lut[0][0]);
            state_arrs[bx](iv, c_h12) = hij_x_arrs[bx](i, j, k, lut[0][1]);
            state_arrs[bx](iv, c_h13) = hij_x_arrs[bx](i, j, k, lut[0][2]);
            state_arrs[bx](iv, c_h22) = hij_x_arrs[bx](i, j, k, lut[1][1]);
            state_arrs[bx](iv, c_h23) = hij_x_arrs[bx](i, j, k, lut[1][2]);
            state_arrs[bx](iv, c_h33) = hij_x_arrs[bx](i, j, k, lut[2][2]);

            state_arrs[bx](iv, c_A11) = Aij_x_arrs[bx](i, j, k, lut[0][0]);
            state_arrs[bx](iv, c_A12) = Aij_x_arrs[bx](i, j, k, lut[0][1]);
            state_arrs[bx](iv, c_A13) = Aij_x_arrs[bx](i, j, k, lut[0][2]);
            state_arrs[bx](iv, c_A22) = Aij_x_arrs[bx](i, j, k, lut[1][1]);
            state_arrs[bx](iv, c_A23) = Aij_x_arrs[bx](i, j, k, lut[1][2]);
            state_arrs[bx](iv, c_A33) = Aij_x_arrs[bx](i, j, k, lut[2][2]);
        }
    });
}

#endif /* INFLATIONEXTRACTION_IMPL_HPP_ */