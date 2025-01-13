/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */


#if !defined(RANDOMFIELD_HPP_)
#error "This file should only be included via RandomField.hpp"
#endif

#ifndef RANDOMFIELD_IMPL_HPP_
#define RANDOMFIELD_IMPL_HPP_

// Used for symmetry condition
inline int RandomField::flip_index(int indx) { return std::abs(N - indx); }

// Used for symmetry condition and calculation of kmag
inline int RandomField::invert_index(int indx) { return (int)(N/2 - std::abs(N/2 - indx)); }

// Used for calculation of polarisation tensors
inline int RandomField::invert_index_with_sign(int indx) 
{ 
    if(indx <= N/2) { return indx; }
    else { return std::abs(N/2 - indx) - N/2; }
}

/****
    Initialisation routines
****/

inline GpuComplex<Real> RandomField::calculate_mode_function(double km, std::string spec_type)
{
    // Deals with k=0 case, undefined if m=0
    if(km < 1.e-23) { return 0.; }
    
    // Stores modulus and argument 
    Real ms_mag = 0.;
    Real ms_arg = 0.;

    // Hubble at t=0, needed for tensor solution
    Real H0 = sqrt((4.0 * M_PI/3.0/pow(m_params.Mp, 2.))
                * (pow(m_background_params.m * m_background_params.phi0, 2.0) 
                    + pow(m_background_params.Pi0, 2.)));

    double kpr = km/H0;
    if (spec_type == "position") // Position mode funcion
    {
        ms_mag = sqrt((1.0/km + H0*H0/pow(km, 3.))/2.);
        ms_arg = atan2((cos(kpr) + kpr*sin(kpr)), (kpr*cos(kpr) - sin(kpr)));
    }
    else if (spec_type == "velocity") // Velocity mode funcion
    {
        ms_mag = sqrt(km/2.);
        ms_arg = -atan2(cos(kpr), sin(kpr));
    }
    else { Error("RandomField::calculate_mode_function Value of spec_type not allowed."); }

    // Construct the mode function and return it
    GpuComplex<Real> ps(ms_mag * cos(ms_arg), ms_mag * sin(ms_arg));
    return ps;
}

inline GpuComplex<Real> RandomField::calculate_random_field(int i, int J, int K, std::string spectrum_type, 
                                                                Real rand_amp, Real rand_phase)
{
    // Storage for the returned value
    GpuComplex<Real> value(0., 0.);

    // Find kmag with FFTW-style inversion on the first two indices
    int j = invert_index(J);
    int k = invert_index(K);
    double kmag = std::sqrt(i*i + j*j + k*k) * 2 * M_PI / m_params.L;

    // Find the linearised solution
    value = calculate_mode_function(kmag, spectrum_type);

    // Add stochastic perturbations
    if(m_params.use_rand == 1)
    {
        BL_PROFILE("RandomField::calculate_random_field Random initialisation is used");

        // Make one random draw for the amplitude and phase
        Real rand_mod = sqrt(-2. * log(rand_amp)); // Rayleigh distribution about |h|
        Real rand_arg = 2. * M_PI * rand_phase;      // Uniform random phase

        // Multiply amplitude by Rayleigh draw
        value *= rand_mod;

        // Apply the random phase, assuming MS phase is accounted for
        Real new_real = value.real() * cos(rand_arg) - value.imag() * sin(rand_arg);
        Real new_imag = value.real() * sin(rand_arg) + value.imag() * cos(rand_arg);
        GpuComplex<Real> new_value(new_real, new_imag);
	
        value = new_value;
    }

    // Apply a window function
    if(m_params.use_window == 1) 
    { 
        BL_PROFILE("RandomField::calculate_random_field Window function is used")
        double ks = m_params.kstar * 2. * M_PI/m_params.L;
        double Dt = m_params.L/m_params.Delta;
        value *= 0.5 * (1. - tanh(Dt * (kmag - ks))); 
    }

    return value;
}

inline Vector<Real> RandomField::calculate_basis_vector(int i, int J, int K, int which_vector)
{
    // Find kmag with FFTW-style inversion on the first two indices
    int j = invert_index_with_sign(J);
    int k = invert_index_with_sign(K);

    Vector<Real> mhat(3, 0.);
    Vector<Real> nhat(3, 0.);

    if (i > 0.) 
    {
        if (k == 0. && j == 0.) { mhat[0] = 1.; mhat[1] = 0.; mhat[2] = 0.; 
                                  nhat[0] = 0.; nhat[1] = 1.; nhat[2] = 0.; 
                                }

        else { mhat[0] = j/sqrt(k*k+j*j); mhat[1] = -k/sqrt(k*k+j*j); mhat[2] = 0.L;
               nhat[0] = k*i/sqrt(i*i*(k*k + j*j) + pow(k*k + j*j, 2.));
               nhat[1] = i*j/sqrt(i*i*(k*k + j*j) + pow(k*k + j*j, 2.));
               nhat[2] = -(k*k + j*j)/sqrt(i*i*(k*k + j*j) + pow(k*k + j*j, 2.)); 
             }
    }

    else if (std::abs(j) > 0) { mhat[0] = 0.; mhat[1] = 0.; mhat[2] = -1.;
                      nhat[0] = -j/sqrt(j*j + k*k);
                      nhat[1] = k/sqrt(j*j + k*k);
                      nhat[2] = 0.; 
                    }

    else if (std::abs(k) > 0) { mhat[0] = 0.; mhat[1] = 1.; mhat[2] = 0.;
                      nhat[0] = 0.; nhat[1] = 0.; nhat[2] = 1.;
                    }

    else if (i==0 && j==0 && k==0) { ; }

    else 
    {
        Error("RandomField::calculate_polarisation_tensors Part of Fourier grid not covered.");
    }

    if(which_vector == 0) { return mhat; }
    else if(which_vector == 1) { return nhat; }
    else { Error("RandomField::calculate_basis_vector Incompatable vector type."); }
}

inline GpuComplex<Real> RandomField::calculate_tensor_initial_conditions(int i, int J, int K, int l, int p, 
                                        GpuComplex<Real> plus_field, GpuComplex<Real> cross_field)
{
    Vector<Real> mhat(3, 0.);
    Vector<Real> nhat(3, 0.);

    mhat = calculate_basis_vector(i, J, K, 0);
    nhat = calculate_basis_vector(i, J, K, 1);

    Real eplus = mhat[l]*mhat[p] - nhat[l]*nhat[p];
    Real ecross = mhat[l]*nhat[p] + nhat[l]*mhat[p];

    return (eplus * plus_field + ecross * cross_field)/std::sqrt(2.);
}

inline void RandomField::apply_nyquist_conditions(int i, int j, int k, Array4<GpuComplex<Real>> const& field)
{
    // Nyquist node condition
    if ((i==0 || i==N/2) && (j==0 || j==N/2) && (k==0 || k== N/2))
    {
        for(int comp = 0; comp < field.nComp(); comp++)
        {
            GpuComplex<Real> temp(field(i, j, k, comp).real(), 0.);
            field(i, j, k, comp) = temp;
        }
    }

    // Nyquist axis condition
    if (i==0 || i==N/2) 
    {
        if((k>N/2 && j==N/2) || (k==0 && j>N/2) || (k>N/2 && j==0) || (k==N/2 && j>N/2))
        {
            for(int comp = 0; comp < field.nComp(); comp++) 
            {
                GpuComplex<Real> temp(field(i, invert_index(j), invert_index(k), comp).real(), 
                                        -field(i, invert_index(j), invert_index(k), comp).imag());
                field(i, j, k, comp) = temp;
            }
        }
        else if(j > N/2)
        {
            for(int comp = 0; comp < field.nComp(); comp++) 
            {
                GpuComplex<Real> temp(field(i, invert_index(j), flip_index(k), comp).real(), 
                                        -field(i, invert_index(j), flip_index(k), comp).imag());
                field(i, j, k, comp) = temp;
            }
        }
    }
}

inline bool RandomField::is_ghost_index(IntVect vector)
{
    bool ret = false;
    for(int d=0; d<3; d++) 
    { 
        if(vector[0] < 0 || vector[0] > N-1) { ret = true; }
    }
    return ret;
}

inline void RandomField::init(amrex::MultiFab &state)
{
    BL_PROFILE("RandomField::init");
    InitRandom(m_params.random_seed);

    BoxArray sba = state.boxArray();
    DistributionMapping sdm = state.DistributionMap();

    // Set up the problem domain and MF ingredients (Real space)
    IntVect domain_low(0, 0, 0);
    IntVect domain_high(N-1, N-1, N-1);
    Box domain(domain_low, domain_high);

    // Make the fft and store the problem domain and MF ingredients (Fourier space)
    FFT::R2C<Real> random_field_fft(domain);
    auto const& [kba, kdm] = random_field_fft.getSpectralDataLayout();

    // Set up the arrays to store the in/out data sets
    cMultiFab hs_k(kba, kdm, 2, 0);
    cMultiFab As_k(kba, kdm, 2, 0);
    cMultiFab hij_k(kba, kdm, 6, 0);
    cMultiFab Aij_k(kba, kdm, 6, 0);

    MultiFab hij_x(sba, sdm, 6, 0);
    MultiFab Aij_x(sba, sdm, 6, 0);

    std::string Filename = "/nfs/st01/hpc-gr-epss/eaf49/GRTeclyn-dump/GRTeclyn-hij-k";
    
    // Loop to create Fourier-space tensor object
    for (MFIter mfi(hs_k); mfi.isValid(); ++mfi) 
    {
        // Make a pointer to the mode functions at this MF box
        Array4<GpuComplex<Real>> const& hs_ptr = hs_k.array(mfi);
        Array4<GpuComplex<Real>> const& hij_ptr = hij_k.array(mfi);

        Array4<GpuComplex<Real>> const& As_ptr = As_k.array(mfi);
        Array4<GpuComplex<Real>> const& Aij_ptr = Aij_k.array(mfi);

        const Box& bx = mfi.fabbox();

        // Loop to create mode functions then hij(k)
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // Find the mode function realisation
            for(int p=0; p<2; p++)
            {
                Real draw1 = amrex::Random();
                Real draw2 = amrex::Random();

                hs_ptr(i, j, k, p) = calculate_random_field(i, j, k, "position", draw1, draw2);
                As_ptr(i, j, k, p) = calculate_random_field(i, j, k, "velocity", draw1, draw2);
            }

            // Find basis tensors and initial tensor realisation
            for (int l=0; l<3; l++) for (int p=0; p<3; p++)
            {
                hij_ptr(i, j, k, lut[l][p]) = calculate_tensor_initial_conditions(i, j, k, l, p, 
                                                hs_ptr(i, j, k, 0), hs_ptr(i, j, k, 1));
                Aij_ptr(i, j, k, lut[l][p]) = calculate_tensor_initial_conditions(i, j, k, l, p, 
                                                As_ptr(i, j, k, 0), As_ptr(i, j, k, 1));
            }
        });

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            //apply_nyquist_conditions(i, j, k, hs_ptr);
            apply_nyquist_conditions(i, j, k, hij_ptr);
            apply_nyquist_conditions(i, j, k, Aij_ptr);

            /*IntVect iv{i, j, k};
            if(i == 0 && j == 0 && k == 1)
            {
                std::cout << "In generation\n";
                std::cout << iv << ": ";
                std::cout << "Fields: " << hs_ptr(i, j, k, 0) << "," << hs_ptr(i, j, k, 1) << "\n";
                Error();
            }*/
        });
    }

    for(int fcomp = 0; fcomp < hij_k.nComp(); fcomp++)
    {
        cMultiFab hij_k_slice(hij_k, make_alias, fcomp, 1);
        MultiFab hij_x_slice(hij_x, make_alias, fcomp, 1);
        random_field_fft.backward(hij_k_slice, hij_x_slice);

        cMultiFab Aij_k_slice(Aij_k, make_alias, fcomp, 1);
        MultiFab Aij_x_slice(Aij_x, make_alias, fcomp, 1);
        random_field_fft.backward(Aij_k_slice, Aij_x_slice);
    }

    /*std::cout << "In generation\n";
    std::cout << "Max: " << hij_x.max(0) << ": ";
    std::cout << "Min: " << hij_x.min(0) << "\n";
    std::cout << "Norm: " << norm << "\n";
    std::cout << "--------\n";*/
    //Error();

    hij_x.mult(norm);
    Aij_x.mult(norm);

    for (int l=0; l<3; l++) { hij_x.plus(1., lut[l][l], 1); }
    Aij_x.mult(-0.5);

    for (MFIter mfi(hij_x); mfi.isValid(); ++mfi) 
    {
        Array4<Real> const& state_ptr = state.array(mfi);
        Array4<Real> const& hij_ptr = hij_x.array(mfi);
        Array4<Real> const& Aij_ptr = Aij_x.array(mfi);

        const Box& bx = mfi.fabbox();
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const IntVect iv{i, j, k};
            bool in_ghost_index = is_ghost_index(iv);
            if(!in_ghost_index)
            {
                state_ptr(iv, c_h11) = hij_ptr(i, j, k, lut[0][0]);
                state_ptr(iv, c_h12) = hij_ptr(i, j, k, lut[0][1]);
                state_ptr(iv, c_h13) = hij_ptr(i, j, k, lut[0][2]);
                state_ptr(iv, c_h22) = hij_ptr(i, j, k, lut[1][1]);
                state_ptr(iv, c_h23) = hij_ptr(i, j, k, lut[1][2]);
                state_ptr(iv, c_h33) = hij_ptr(i, j, k, lut[2][2]);

                state_ptr(iv, c_A11) = Aij_ptr(i, j, k, lut[0][0]);
                state_ptr(iv, c_A12) = Aij_ptr(i, j, k, lut[0][1]);
                state_ptr(iv, c_A13) = Aij_ptr(i, j, k, lut[0][2]);
                state_ptr(iv, c_A22) = Aij_ptr(i, j, k, lut[1][1]);
                state_ptr(iv, c_A23) = Aij_ptr(i, j, k, lut[1][2]);
                state_ptr(iv, c_A33) = Aij_ptr(i, j, k, lut[2][2]);
            }
        });
    }
}

/****
    Extraction routines
****/

inline void RandomField::print_tensor_moment(int moment_order, MultiFab &field)
{
    std::cout << field.nComp() << "\n";
    Error();
}

inline void RandomField::print_power_spectrum(cMultiFab &field_array, SmallDataIO &power_spec_file, int component)
{ 
    double kiso_max = std::sqrt(3.) * N * M_PI / m_params.L; // maximum unitful k vector length on the isotropic (diagonal) axis
    double dkiso = sqrt(3.)*2.*M_PI/m_params.L;          // step in k across the isotropic axis
    double tolerance = 1.e-12;

    // check the stepping along the diagonal is consistent
    if (kiso_max/dkiso - N/2 > tolerance)
    {
        Error("RandomField::print_power_spectrum Isotropic k axis is too large.");
    }
    // check you aren't sampling above the max sampling rate
    else if (m_params.bin_number > kiso_max/dkiso)
    {
        Error("RandomField::print_power_spectrum Bin number is too large.");
    }

    // Set up isotropic k axis
    double dk_to_bin = (double)m_params.bin_number/((double)N/2);
    double kmag = 0.;
    Vector<Real> kiso(N/2+1, 0.);
    for (int s=0; s<=N/2; s++) { kiso[s] = s*dkiso; }

    Vector<Real> power_spec(m_params.bin_number + 1, 0.);
    Vector<Real> k_counter(m_params.bin_number + 1, 0.);

    // Loop to create Fourier-space tensor object
    MFIter::allowMultipleMFIters(true);
    for (MFIter mfi(field_array); mfi.isValid(); ++mfi) 
    {
        // Make a pointer to the mode functions at this MF box
        Array4<GpuComplex<Real>> const& field_ptr = field_array.array(mfi);
        const Box& bx = mfi.fabbox();

        // Loop to create mode functions then hij(k)
        amrex::ParallelFor(bx, [=, &power_spec, &k_counter] AMREX_GPU_DEVICE (int i, int J, int K) noexcept
        {
            int j = invert_index(J);
            int k = invert_index(K);
            IntVect iv{i, j, k};

            double kmag = std::sqrt(i*i + j*j + k*k) * 2 * M_PI / m_params.L;

            // make sure you're still in the domain
            if(kmag > kiso_max) 
            { 
                std::cout << iv << "\n";
                std::cout << kmag << "," << kiso_max << "\n";
                Error("RandomField::print_power_spectrum Found magnitude larger than (N/2,N/2,N/2)."); 
            }

            for (int s=1; s<=N/2; s++) 
            {
                // If smaller than the smallest bin
                if(kmag < kiso[0])
                {
                    std::cout << iv << "\n";
                    Error("RandomField::print_power_spectrum kmag below the kiso domain.");
                }

                // If you're larger than the largest bin
                else if(kmag - kiso[N/2] > tolerance)
                {
                    std::cout << iv << "\n";
                    Error("RandomField::print_power_spectrum kmag above the kiso domain.");
                }

                // If you're somewhere in the middle
                else if (kmag < kiso[s] && kmag >= kiso[(s-1)]) 
                {
                    int bin_index = round((s-1)*dk_to_bin);
                    if(bin_index > m_params.bin_number)
                    {
                        std::cout << iv << "\n";
                        Error("RandomField::print_power_spectrum Bin index larger than requested bin number.");
                    }

                    Real power = (std::pow(field_ptr(i, j, k, component).real(), 2.0) 
                                            + std::pow(field_ptr(i, j, k, component).imag(), 2.0));

                    power_spec[bin_index] += power;
                    k_counter[bin_index] += 1;

                    break;
                }

                else if(s == N/2)
                { 
                    Real power = (std::pow(field_ptr(i, j, k, component).real(), 2.0) 
                                            + std::pow(field_ptr(i, j, k, component).imag(), 2.0));

                    power_spec[m_params.bin_number + 1] += power;
                    k_counter[m_params.bin_number + 1] += 1;

                    break;
                }

                // If you've reached the largest bin but not been captured
                else if(s > N/2)
                { 
                    std::cout << iv << "\n";
                    std::cout << kmag << "\n";
                    std::cout << kiso[s] << "," << kiso[s-1] << "\n";
                    Error("RandomField::print_power_spectrum Part of the spectrum isn't captured.");
                }

                // If you haven't found the right bin yet
                else { continue; }
            }
        });
    }

    for(int s=0; s<=N/2; s++)
    {
        if(k_counter[s] != 0) { power_spec_file.write_data_line({kiso[s], power_spec[s]/k_counter[s]}); }
    }
}

inline void RandomField::extract(MultiFab &state, std::string data_path, Real dt, Real cur_time, int restart_time, int first_step)
{
    BL_PROFILE("RandomField::init");

    // Extract MultiFab ingredients from state
    BoxArray sba = state.boxArray();
    DistributionMapping sdm = state.DistributionMap();
    MultiFab hij_x(sba, sdm, 6, 0);

    Copy(hij_x, state, c_h11, lut[0][0], 1, 0);
    Copy(hij_x, state, c_h12, lut[0][1], 1, 0);
    Copy(hij_x, state, c_h13, lut[0][2], 1, 0);
    Copy(hij_x, state, c_h22, lut[1][1], 1, 0);
    Copy(hij_x, state, c_h23, lut[1][2], 1, 0);
    Copy(hij_x, state, c_h33, lut[2][2], 1, 0);

    // Undo the normalisation and delta function 
    for (int l=0; l<3; l++) { hij_x.plus(-1., lut[l][l], 1); }
    hij_x.mult(1./norm);

    /*std::cout << "In extraction\n";
    std::cout << "Max: " << hij_x.max(0) << ": ";
    std::cout << "Min: " << hij_x.min(0) << "\n";
    std::cout << "Norm: " << norm << "\n";
    std::cout << "--------\n";
    Error();*/

    /*if(i != 0 || j != 0 || k != 0)
    {
        std::cout << iv << ": ";
        std::cout << "Fields: " << hij_ptr(i, j, k, 0) << "," << hij_ptr(i, j, k, 1) << "\n";
        Error();
    }*/

    // Set up the problem domain and MF ingredients (Fourier space)
    IntVect domain_low(0, 0, 0);
    IntVect domain_high(N-1, N-1, N-1);
    Box domain(domain_low, domain_high);
    FFT::R2C<Real> random_field_fft(domain);
    auto const& [kba, kdm] = random_field_fft.getSpectralDataLayout();

    // Set up the arrays to store the Fourier data sets
    cMultiFab hs_k(kba, kdm, 2, 0);
    cMultiFab hij_k(kba, kdm, 6, 0);

    for(int fcomp = 0; fcomp < hij_x.nComp(); fcomp++)
    {
        cMultiFab hij_k_slice(hij_k, make_alias, fcomp, 1);
        MultiFab hij_x_slice(hij_x, make_alias, fcomp, 1);
        random_field_fft.forward(hij_x_slice, hij_k_slice);
    }

    /*std::cout << "In extraction\n";
    std::cout << "Max: " << hij_k.max(0) << ": ";
    std::cout << "Min: " << hij_k.min(0) << "\n";
    std::cout << "Norm: " << norm << "\n";
    std::cout << "--------\n";
    Error();*/

    for(int comp = 0; comp < 6; comp++)
    {
        hij_k.mult(std::pow(N, -3.), comp, 1); 
    }

    // Loop to create Fourier-space tensor object
    for (MFIter mfi(hij_k); mfi.isValid(); ++mfi) 
    {
        // Make a pointer to the mode functions at this MF box
        Array4<GpuComplex<Real>> const& hs_ptr = hs_k.array(mfi);
        Array4<GpuComplex<Real>> const& hij_ptr = hij_k.array(mfi);
        const Box& bx = mfi.fabbox();

        // Loop to create mode functions then hij(k)
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            IntVect iv{i, j, k};
            Vector<Real> mhat(3, 0.);
            Vector<Real> nhat(3, 0.);

            mhat = calculate_basis_vector(i, j, k, 0);
            nhat = calculate_basis_vector(i, j, k, 1);

            Real eplus = 0.;
            Real ecross = 0.;

            // Find basis tensors and initial tensor realisation
            for (int l=0; l<3; l++) for (int p=0; p<3; p++)
            {
                //hij_k(i, j, k, lut[l][p]) *= 1./std::pow(N, 3.);

                eplus = mhat[l]*mhat[p] - nhat[l]*nhat[p];
                ecross = mhat[l]*nhat[p] + nhat[l]*mhat[p];

                hs_ptr(i, j, k, 0) += (hij_ptr(i, j, k, lut[l][p]) * eplus)/std::sqrt(2.);
                hs_ptr(i, j, k, 1) += (hij_ptr(i, j, k, lut[l][p]) * ecross)/std::sqrt(2.);
            }

            /*if(i == 0 && j == 0 && k == 1)
            {
                std::cout << "In extraction\n";
                std::cout << iv << ": ";
                std::cout << "Fields: " << hs_ptr(i, j, k, 0) << "," << hs_ptr(i, j, k, 1) << "\n";
                Error();
            }*/
            
        });

        /*amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            //apply_nyquist_conditions(i, j, k, hs_ptr);
            apply_nyquist_conditions(i, j, k, hij_ptr);
            apply_nyquist_conditions(i, j, k, Aij_ptr);
        });*/

        // THERE MAY BE AN MPI ISSUE HERE
        /*if(m_params.calc_binned_power_spectrum) 
        {
            for(int comp = 0; comp < hs_k.nComp(); comp++)
            {
                SmallDataIO spectrum_file(data_path+"spectrum-comp-"+std::to_string(comp)+"-time-", 
                                            dt, cur_time, restart_time, SmallDataIO::NEW, first_step, ".dat");
                print_power_spectrum(hs_k, spectrum_file, comp);
            }
        }*/
    }
}

#endif /* RANDOMFIELD_IMPL_HPP_*/
