/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */


#if !defined(RANDOMFIELDINIT_HPP_)
#error "This file should only be included via RandomFieldInit.hpp"
#endif

#ifndef RANDOMFIELDINIT_IMPL_HPP_
#define RANDOMFIELDINIT_IMPL_HPP_

// Generate unique random draws for each MFI box.
inline void RandomFieldInit::make_random_draws(MultiFab &rand_fab, const Box &domain, 
                                               const int seed)
{
    BoxArray ba = rand_fab.boxArray();
    DistributionMapping dm = rand_fab.DistributionMap();
    MultiFab tmp(ba, dm, rand_fab.nComp(), 0, MFInfo{}.SetArena(The_Cpu_Arena()));

    // CHECKME: should this go inside the MFI loop, or does it not matter?
    std::mt19937 generator(seed);
    std::uniform_real_distribution<Real> distribution(Real(0), Real(1));

    for(MFIter mfi(tmp); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.validbox();
        auto const& tmp_ptr = tmp.array(mfi);

        auto offset = domain.index(bx.smallEnd()) * tmp.nComp();
        for(int ofs = 0; ofs < offset; ofs++)
        {
            distribution(generator);
        }
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k)
        {
            for(int l=0; l<tmp.nComp(); l++)
            {
                tmp_ptr(i, j, k, l) = distribution(generator);
            }
        });
    }

    rand_fab.ParallelCopy(tmp);
}

inline GpuComplex<Real> RandomFieldInit::find_in_stoiic(const double km, const int field_indx, 
                                                    const std::string field_type)
{
    // Assume no average
    if(km == 0) { return GpuComplex<Real>{0., 0.}; }

    // Find the index where this k appears
    int spec_index;
    for(int idx = 0; idx < m_params.init_k.size(); idx++)
    {
        if(std::abs(km - m_params.init_k[idx]) < 1e-13) { spec_index = idx; break; }
        else if (idx == m_params.init_k.size() - 1) 
        { 
            Print() << km << "\n"; 
            Error("The above k was not found in the STOIIC file."); 
        }
    }

    // Return the field at this k
    if(field_type == "tensor")
    {
        return (GpuComplex<Real>{m_params.tensor_ps[2*field_indx][spec_index], 
                                 m_params.tensor_ps[2*field_indx+1][spec_index]});
    }
    else if(field_type == "scalar")
    {
        return (GpuComplex<Real>{m_params.scalar_ps[2*field_indx][spec_index], 
                                 m_params.scalar_ps[2*field_indx+1][spec_index]});
    }
    else 
    { 
        Error("RandomFieldInit::find_in_stoiic field cannot be found."); 
        return GpuComplex<Real>{0., 0.}; 
    }
}

// Returns analytic power spectrum in modulus/argument form
inline GpuComplex<Real> RandomFieldInit::calculate_mode_function(const double km, 
                                                             const int spec_indx)
{
    // Deals with k=0 case, which is undefined if m=0
    if(km < 1.e-15) { return GpuComplex<Real>{0., 0.}; }
    
    // Stores modulus and argument 
    Real ms_mag = 0.;
    Real ms_arg = 0.;

    double kpr = km/H0;
    if (spec_indx == 0) // Position mode funcion
    {
        ms_mag = sqrt((1.0/km + H0*H0/pow(km, 3.))/2.);
        ms_arg = atan2((cos(kpr) + kpr*sin(kpr)), (kpr*cos(kpr) - sin(kpr)));
    }
    else if (spec_indx == 1) // Velocity mode funcion
    {
        ms_mag = sqrt(km/2.);
        ms_arg = -atan2(cos(kpr), sin(kpr));
    }
    else { Error("RandomFieldInit::calculate_mode_function Value of spec_type not allowed."); }

    // Construct the mode function and return it
    GpuComplex<Real> ps(ms_mag * cos(ms_arg), ms_mag * sin(ms_arg));
    return ps;
}

// Turns analytic PS into GRF and applies window function if requested
inline GpuComplex<Real> RandomFieldInit::calculate_random_field(const IntVect iv, const int field_index, 
                                                            const Real rand_amp, const Real rand_phase, 
                                                            std::string field_type)
{
    GpuComplex<Real> value(0., 0.);

    // Find kmag with FFTW-style inversion on the last two indices
    double kmag = m_params.get_kmag(iv);

    // Find the analytic power spectrum
    if(m_params.read_from_stoiic) 
    { 
        value = find_in_stoiic(kmag, field_index, field_type); 
    }
    else { value = calculate_mode_function(kmag, field_index); }

    // Add stochastic perturbations
    if(m_params.use_rand == 1)
    {
        BL_PROFILE("RandomFieldInit::calculate_random_field Random initialisation is used");

        // Make one random draw for the amplitude and phase individually
        Real rand_mod = sqrt(-2. * log(rand_amp)); // Rayleigh distribution about |h|
        Real rand_arg = 2. * M_PI * rand_phase;      // Uniform random phase

        // Multiply amplitude by Rayleigh draw
        value *= rand_mod;

        // Apply the random phase (assuming MS phase is accounted for)
        Real new_real = value.real() * cos(rand_arg) - value.imag() * sin(rand_arg);
        Real new_imag = value.real() * sin(rand_arg) + value.imag() * cos(rand_arg);
        GpuComplex<Real> new_value(new_real, new_imag);
	
        value = new_value;
    }

    // Apply a window function if requested
    if(m_params.use_window == 1) 
    { 
        BL_PROFILE("RandomFieldInit::calculate_random_field Window function is used")
        double ks = std::sqrt(3.) * m_params.N * M_PI / m_params.L / 5. / 2.;
        //m_params.kstar * 2. * M_PI/m_params.L;
        double Dt = m_params.L/m_params.Delta;
        value *= 0.5 * (1.0 - tanh(Dt * (kmag - ks))); 
    }

    return value;
}

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
    IntVect k_domain_high(m_params.N/2, m_params.N-1, m_params.N-1);
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
    IntVect x_domain_high(m_params.N-1, m_params.N-1, m_params.N-1);
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
    const auto &hij_k_arrs = hij_k.arrays();
    const auto &As_arrs = As_k.arrays();
    const auto &Aij_k_arrs = Aij_k.arrays();
    const auto &scalar_field_arrs = scalar_fields_k.arrays();

    Print() << "Starting initial condition generation/read in...\n";
    
    ParallelFor(hs_k,
        [=] AMREX_GPU_DEVICE (int bx, int i, int j, int k)
    {
        IntVect iv = {i, j, k};
        if(iv != IntVect{0, 0, 0})
        {
            if(m_params.scalar_init)
            {
                for(int f=0; f<4; f++)
                {
                    Real draw1 = scalar_draw_arrs[bx](i, j, k, 0);
                    Real draw2 = scalar_draw_arrs[bx](i, j, k, 1);

                    scalar_field_arrs[bx](i, j, k, f) = calculate_random_field(iv, f, draw1, draw2, "scalar");
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
                
                // Construct polarisation tensors from basis vectors
                Tensor<2, Real> eplus = m_params.calculate_polarisation_tensor(iv, 0);
                Tensor<2, Real> ecross = m_params.calculate_polarisation_tensor(iv, 1);

                // Find basis tensors and initial tensor realisation
                for (int l=0; l<3; l++) for (int p=0; p<3; p++)
                {
                    hij_k_arrs[bx](i, j, k, InflationUtils::lut[l][p]) = (hs_arrs[bx](i, j, k, 0) * eplus[l][p] 
                                                                        + hs_arrs[bx](i, j, k, 1) * ecross[l][p]);
                    Aij_k_arrs[bx](i, j, k, InflationUtils::lut[l][p]) = (As_arrs[bx](i, j, k, 0) * eplus[l][p] 
                                                                        + As_arrs[bx](i, j, k, 1) * ecross[l][p]);
                }
            }
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
    Aij_x.mult(-0.5);

    const auto &state_arrs = state.arrays();
    const auto &hij_x_arrs = hij_x.arrays();
    const auto &Aij_x_arrs = Aij_x.arrays();
    const auto &scalar_field_x_arrs = scalar_fields_x.arrays();

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
            state_arrs[bx](iv, c_h11) += hij_x_arrs[bx](i, j, k, InflationUtils::lut[0][0]);
            state_arrs[bx](iv, c_h12) += hij_x_arrs[bx](i, j, k, InflationUtils::lut[0][1]);
            state_arrs[bx](iv, c_h13) += hij_x_arrs[bx](i, j, k, InflationUtils::lut[0][2]);
            state_arrs[bx](iv, c_h22) += hij_x_arrs[bx](i, j, k, InflationUtils::lut[1][1]);
            state_arrs[bx](iv, c_h23) += hij_x_arrs[bx](i, j, k, InflationUtils::lut[1][2]);
            state_arrs[bx](iv, c_h33) += hij_x_arrs[bx](i, j, k, InflationUtils::lut[2][2]);

            state_arrs[bx](iv, c_A11) += Aij_x_arrs[bx](i, j, k, InflationUtils::lut[0][0]);
            state_arrs[bx](iv, c_A12) += Aij_x_arrs[bx](i, j, k, InflationUtils::lut[0][1]);
            state_arrs[bx](iv, c_A13) += Aij_x_arrs[bx](i, j, k, InflationUtils::lut[0][2]);
            state_arrs[bx](iv, c_A22) += Aij_x_arrs[bx](i, j, k, InflationUtils::lut[1][1]);
            state_arrs[bx](iv, c_A23) += Aij_x_arrs[bx](i, j, k, InflationUtils::lut[1][2]);
            state_arrs[bx](iv, c_A33) += Aij_x_arrs[bx](i, j, k, InflationUtils::lut[2][2]);
        }
    });
}

#endif /* RANDOMFIELDINIT_IMPL_HPP_ */