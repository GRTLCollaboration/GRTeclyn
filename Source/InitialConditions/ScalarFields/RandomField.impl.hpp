/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */


#if !defined(RANDOMFIELD_HPP_)
#error "This file should only be included via RandomField.hpp"
#endif

#ifndef RANDOMFIELD_IMPL_HPP_
#define RANDOMFIELD_IMPL_HPP_

/****
    Small functions (less than 10 lines)
****/

// Nyquist condition
inline int RandomField::flip_index(const int indx) { return std::abs(N - indx); }

// Nyquist condition and calculation of kmag
inline int RandomField::invert_index(const int indx) { return (int)(N/2 - std::abs(N/2 - indx)); }

// For calculation of polarisation tensors
inline int RandomField::invert_index_with_sign(const int indx) 
{ 
    if(indx <= N/2) { return indx; }
    else { return std::abs(N/2 - indx) - N/2; }
}

inline Real RandomField::get_kmag(IntVect iv)
{
    const int i = iv[0];
    const int j = invert_index(iv[1]);
    const int k = invert_index(iv[2]);
    return std::sqrt(i*i + j*j + k*k) * 2 * M_PI / m_params.L;
}

// Ensures no calculation on ghost cells
inline bool RandomField::is_ghost_index(const IntVect vector)
{
    bool ret = false;
    for(int d=0; d<3; d++) 
    { 
        if(vector[d] < 0 || vector[d] > N-1) { ret = true; }
    }
    return ret;
}

// Makes subdirectories in data/
inline std::string RandomField::make_subdirectory(const std::string base, const std::string dir, const int is_first_step)
{
    std::string new_path = base+"../"+dir+"/";
    if(is_first_step)
    {
        if (FilesystemTools::directory_exists(base)) { FilesystemTools::mkdir_recursive(new_path); }
        else 
        { 
            Print() << "Directory creation failed for " << new_path << "\n";
            Error("RandomField::extract Data directory has not been created."); 
        }
    }
    return new_path;
}

// Creates a custom data file layout 
inline void RandomField::assign_statistics_data(Vector<std::string> &header_storage, const std::string name, 
                            Vector<Real> &data_storage, const Vector<Real> data, const int component, const int num_comps,
                            const Vector<int>::const_iterator itr, const Vector<int>::const_iterator start, 
                            const int is_first_step)
{
    int loc = component + num_comps*(itr - start);
    if(is_first_step) 
    { 
        header_storage[loc] =  name; 
    }
    data_storage[loc] = data[component];
}

/****
    Tests
****/

inline void RandomField::Test_is_trace_free(MultiFab &field)
{
    for (MFIter mfi(field); mfi.isValid(); ++mfi) 
    {
        Array4<Real> const& field_ptr = field.array(mfi);
        const Box& bx = mfi.fabbox();

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            IntVect iv{i, j, k};
            Real sum = 0.;

            for(int l=0; l<3; l++)
            {
                sum += field_ptr(i, j, k, lut[l][l]);
            }

            if(std::abs(sum) > tolerance)
            {
                Print() << iv << ": " << sum;
                Error("RandomField::Test_is_trace_free Trace-free test failed here.");
            }
        });
    }
    
}

/****
    Initialisation routines
****/

// Generate unique random draws for each MFI box.
inline void RandomField::make_random_draws(MultiFab &rand_fab, Box &domain, const int seed)
{
    BoxArray ba = rand_fab.boxArray();
    DistributionMapping dm = rand_fab.DistributionMap();
    MultiFab tmp(ba, dm, 6, 0, MFInfo{}.SetArena(The_Cpu_Arena()));

    for(MFIter mfi(tmp); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.validbox();
        auto const& tmp_ptr = tmp.array(mfi);

        std::mt19937 generator(seed);
        std::uniform_real_distribution<Real> distribution(Real(0), Real(1));

        auto offset = domain.index(bx.smallEnd()) * 6;
        for(int ofs = 0; ofs < offset; ofs++)
        {
            distribution(generator);
        }
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k)
        {
            for(int l=0; l<6; l++)
            {
                tmp_ptr(i, j, k, l) = distribution(generator);
            }
        });
    }

    rand_fab.ParallelCopy(tmp);
}

// Returns analytic power spectrum in modulus/argument form
inline GpuComplex<Real> RandomField::calculate_mode_function(const double km, const int spec_indx)
{
    // Deals with k=0 case, which is undefined if m=0
    if(km < 1.e-23) { return 0.; }
    
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
    else { Error("RandomField::calculate_mode_function Value of spec_type not allowed."); }

    // Construct the mode function and return it
    GpuComplex<Real> ps(ms_mag * cos(ms_arg), ms_mag * sin(ms_arg));
    return ps;
}

inline GpuComplex<Real> RandomField::find_in_stoiic(const double km, const int field_indx, std::string field_type)
{
    GpuComplex<Real> zero(0., 0.);
    if(km == 0) { return zero; }

    int spec_index;
    for(int idx = 0; idx < m_params.init_k.size(); idx++)
    {
        if(std::abs(km - m_params.init_k[idx]) < 1e-15) { spec_index = idx; break; }
        else if (idx == m_params.init_k.size() - 1) { AllPrint() << km << "\n"; Error("The above k was not found in the STOIIC file."); }
    }

    if(field_type == "tensor")
    {
        return GpuComplex<Real>{m_params.tensor_ps[2*field_indx][spec_index], m_params.tensor_ps[2*field_indx+1][spec_index]};
    }
    else if(field_type == "scalar")
    {
        return GpuComplex<Real>{m_params.scalar_ps[2*field_indx][spec_index], m_params.scalar_ps[2*field_indx+1][spec_index]};
    }
    else { Error("RandomField::find_in_stoiic field cannot be found."); return GpuComplex<Real>{0., 0.}; }
}

// Turns analytic PS into GRF and applies window function if requested
inline GpuComplex<Real> RandomField::calculate_random_field(const IntVect iv, const int field_index, 
                                                            const Real rand_amp, const Real rand_phase, 
                                                            std::string field_type)
{
    GpuComplex<Real> value(0., 0.);

    // Find kmag with FFTW-style inversion on the last two indices
    double kmag = get_kmag(iv);

    // Find the analytic power spectrum
    if(m_params.read_from_stoiic) { value = find_in_stoiic(kmag, field_index, field_type); }
    else { value = calculate_mode_function(kmag, field_index); }

    // Add stochastic perturbations
    if(m_params.use_rand == 1)
    {
        BL_PROFILE("RandomField::calculate_random_field Random initialisation is used");

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
        BL_PROFILE("RandomField::calculate_random_field Window function is used")
        double ks = std::sqrt(3.) * N * M_PI / m_params.L / 5. / 2.;//m_params.kstar * 2. * M_PI/m_params.L;
        double Dt = m_params.L/m_params.Delta;
        value *= 0.5 * (1.0 - tanh(Dt * (kmag - ks))); 
    }

    return value;
}

// Calculates basis vectors required for polarisation tensors
inline Vector<Real> RandomField::calculate_basis_vector(const IntVect iv, const int which_vector)
{
    // FFTW-style inversion with sign on the last two indices
    int i = iv[0];
    int j = invert_index_with_sign(iv[1]);
    int k = invert_index_with_sign(iv[2]);

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

    if(m_params.alpha != 0)
    {
        double m_alpha = m_params.alpha * M_PI/180.;
        for (int i=0; i<3; i++)
        {
            mhat[i] = cos(m_alpha) * mhat[i] + sin(m_alpha) * nhat[i];
            nhat[i] = -sin(m_alpha) * mhat[i] + cos(m_alpha) * nhat[i];
        }
    }

    if(which_vector == 0) { return mhat; }
    else if(which_vector == 1) { return nhat; }
    else { Error("RandomField::calculate_basis_vector Incompatable vector type."); return mhat; }
}

// Assembles full tensor initial conditions given two mode functions
inline GpuComplex<Real> RandomField::calculate_tensor_initial_conditions(const IntVect iv, const int l, const int p, 
                                                                         const GpuComplex<Real> plus_field, 
                                                                         const GpuComplex<Real> cross_field)
{
    Vector<Real> mhat(3, 0.);
    Vector<Real> nhat(3, 0.);

    mhat = calculate_basis_vector(iv, 0);
    nhat = calculate_basis_vector(iv, 1);

    // Assemble the polarisation tensors
    Real eplus = mhat[l]*mhat[p] - nhat[l]*nhat[p];
    Real ecross = mhat[l]*nhat[p] + nhat[l]*mhat[p];

    return (eplus * plus_field + ecross * cross_field);
}

// Applies above Nyquist conditions to a given MF
inline void RandomField::apply_nyquist_conditions(cMultiFab &field)
{
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

// Main initialisation routine
inline void RandomField::init(amrex::MultiFab &state)
{
    BL_PROFILE("RandomField::init");

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

    // Make the Fourier transform
    IntVect x_domain_high(N-1, N-1, N-1);
    Box x_domain(domain_low, x_domain_high);
    FFT::R2C<Real> tensor_fft(x_domain, FFT::Info().setBatchSize(hij_k.nComp()));
    FFT::R2C<Real> scalar_fft(x_domain, FFT::Info().setBatchSize(scalar_fields_k.nComp()));

    MultiFab random_draws(kba, kdm, 6, 0);
    make_random_draws(random_draws, k_domain, m_params.random_seed);
    MultiFab tensor_draws(random_draws, amrex::make_alias, 0, 4);
    MultiFab scalar_draws(random_draws, amrex::make_alias, 4, 2);

    Print() << "Starting initial condition generation/read in...\n";

    // std::string Filename = "/cephfs/home/eaf49/GRTeclyn-workspace/Examples/ScalarField/comp-to-dparams.txt";
    for (MFIter mfi(hs_k); mfi.isValid(); ++mfi) 
    {
        // Define the domain on this MPI rank
        const Box& bx = mfi.fabbox();
        auto const& tensor_draw_ptr = tensor_draws.const_array(mfi);
        auto const& scalar_draw_ptr = scalar_draws.const_array(mfi);

        // Make a pointer to the mode functions at this MF box
        Array4<GpuComplex<Real>> const& hs_ptr = hs_k.array(mfi);
        Array4<GpuComplex<Real>> const& hij_ptr = hij_k.array(mfi);

        Array4<GpuComplex<Real>> const& As_ptr = As_k.array(mfi);
        Array4<GpuComplex<Real>> const& Aij_ptr = Aij_k.array(mfi);

        Array4<GpuComplex<Real>> const& scalar_fields_ptr = scalar_fields_k.array(mfi);

        // Loop to create mode functions, then hij(k) and Aij(k)
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            IntVect iv = {i, j, k};

            if(m_params.scalar_init)
            {
                for(int f=0; f<4; f++)
                {
                    Real draw1 = scalar_draw_ptr(i, j, k, 0);
                    Real draw2 = scalar_draw_ptr(i, j, k, 1);

                    scalar_fields_ptr(i, j, k, f) = calculate_random_field(iv, f, draw1, draw2, "scalar");
                }
            }

            if(m_params.tensor_init)
            {
                // Find the mode function realisation
                for(int p=0; p<2; p++)
                {
                    Real draw1 = tensor_draw_ptr(i, j, k, 2*p);
                    Real draw2 = tensor_draw_ptr(i, j, k, 2*p+1);

                    hs_ptr(i, j, k, p) = calculate_random_field(iv, 0, draw1, draw2, "tensor");
                    As_ptr(i, j, k, p) = calculate_random_field(iv, 1, draw1, draw2, "tensor");
                }

                // Find basis tensors and initial tensor realisation
                for (int l=0; l<3; l++) for (int p=0; p<3; p++)
                {
                    hij_ptr(i, j, k, lut[l][p]) = calculate_tensor_initial_conditions(iv, l, p, hs_ptr(i, j, k, 0), hs_ptr(i, j, k, 1));
                    Aij_ptr(i, j, k, lut[l][p]) = calculate_tensor_initial_conditions(iv, l, p, As_ptr(i, j, k, 0), As_ptr(i, j, k, 1));

                    // if(i == 1 && j == 0 && k == 1)
                    // {
                    //     Print() << "In init: \n";
                    //     Print() << l << ", " << p << ": ";
                    //     Print() << hij_ptr(i, j, k, lut[l][p]) << "\n";
                    // }
                }
            }
        });
    }

    // Apply the Nyquist conditions
    apply_nyquist_conditions(hs_k);
    apply_nyquist_conditions(hij_k);
    apply_nyquist_conditions(Aij_k);
    apply_nyquist_conditions(scalar_fields_k);

    // Do the Fourier transform
    tensor_fft.backward(hij_k, hij_x);
    tensor_fft.backward(Aij_k, Aij_x);
    scalar_fft.backward(scalar_fields_k, scalar_fields_x);

    // Apply normalisation into physical units
    hij_x.mult(norm);
    Aij_x.mult(norm);
    scalar_fields_x.mult(norm);

    // Test is trace-free
    Test_is_trace_free(hij_x);
    Test_is_trace_free(Aij_x);

    // Convert to BSSN variables using the BSSN-CPT dictionary
    for (int l=0; l<3; l++) { hij_x.plus(1., lut[l][l], 1); }
    Aij_x.mult(-0.5);

    // Put these initial conditions into state
    for (MFIter mfi(hij_x); mfi.isValid(); ++mfi) 
    {
        Array4<Real> const& state_ptr = state.array(mfi);
        Array4<Real> const& hij_ptr = hij_x.array(mfi);
        Array4<Real> const& Aij_ptr = Aij_x.array(mfi);
        Array4<Real> const& scalar_ptr = scalar_fields_x.array(mfi);

        const Box& bx = mfi.fabbox();
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const IntVect iv{i, j, k};
            bool in_ghost_index = is_ghost_index(iv);
            if(!in_ghost_index)
            {
                if(m_params.tensor_init)
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

                if(m_params.scalar_init)
                {
                    state_ptr(iv, c_phi) += scalar_ptr(i, j, k, 0);
                    state_ptr(iv, c_Pi) += scalar_ptr(i, j, k, 1);
                    state_ptr(iv, c_chi) += scalar_ptr(i, j, k, 2);
                    state_ptr(iv, c_K) += scalar_ptr(i, j, k, 3);
                }
            }
        });
    }
}

/****
    Extraction routines
****/

// Calculates and prints the power spectrum
inline void RandomField::print_power_spectrum(cMultiFab &field_array, SmallDataIO &power_spec_file, const int component = 0)
{ 
    // Set up the isotropic k axis bounds
    double kiso_max = std::sqrt(3.) * N * M_PI / m_params.L;
    double dkiso = sqrt(3.)*2.*M_PI/m_params.L;
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
    // check your bin number isn't greater than the max resolvable bins
    else if(m_params.bin_number > m_params.N_readin/2)
    {
        Error("RandomField::print_power_spectrum bin number must be less than N/2.");
    }

    // Set up isotropic k axis and PS map
    double dk_to_bin = (double)m_params.bin_number/((double)N/2);
    double kmag = 0.;
    Vector<Real> kiso(N/2+1, 0.);

    Vector<Real> ps_map(m_params.bin_number+1, 0.);
    Vector<int> kcount(m_params.bin_number+1, 0);

    for (int s=0; s<=N/2; s++) { kiso[s] = s*dkiso; }

    // Loop to bin the power spectrum at each point
    MFIter::allowMultipleMFIters(true); // Needed to pass the map to the ParallelFor loop
    for (MFIter mfi(field_array); mfi.isValid(); ++mfi) 
    {
        // Make a pointer to the mode function at this MF box
        Array4<GpuComplex<Real>> const& field_ptr = field_array.array(mfi);
        const Box& bx = mfi.fabbox();

        amrex::ParallelFor(bx, [=, &ps_map, &kcount] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // Check to see if you're in a ghost cell
            IntVect iv{i, j, k};
            bool in_ghost_index = is_ghost_index(iv);
            if(!in_ghost_index)
            {
                double kmag = get_kmag(iv);

                // make sure you're still in the domain
                if(kmag - kiso_max > tolerance) 
                { 
                    Print() << iv << "\n";
                    Print() << kmag << "," << kiso_max << "\n";
                    Error("RandomField::print_power_spectrum Found magnitude larger than (N/2,N/2,N/2)."); 
                }

                // Loop over the isotropic axis
                for (int s=1; s<=N/2; s++) 
                {
                    // If smaller than the smallest bin
                    if(kmag < kiso[0])
                    {
                        Print() << iv << "\n";
                        Error("RandomField::print_power_spectrum kmag below the kiso domain.");
                    }

                    // If you're larger than the largest bin
                    else if(kmag - kiso[N/2] > tolerance)
                    {
                        Print() << iv << "\n";
                        Error("RandomField::print_power_spectrum kmag above the kiso domain.");
                    }

                    // If you're somewhere in the middle
                    else if (kmag < kiso[s] && kmag >= kiso[(s-1)]) 
                    {
                        Real power = (std::pow(field_ptr(i, j, k, component).real(), 2.0) 
                                    + std::pow(field_ptr(i, j, k, component).imag(), 2.0));
                        
                        Gpu::Atomic::Add(&kcount[s-1], 1);
                        if(power > tolerance)
                        {
                            Gpu::Atomic::Add(&ps_map[s-1], power);   
                        }

                        break;
                    }

                    // If you're at the largest bin
                    else if (kmag == kiso[N/2])
                    { 
                        Real power = (std::pow(field_ptr(i, j, k, component).real(), 2.0) 
                                    + std::pow(field_ptr(i, j, k, component).imag(), 2.0));
                        
                        Gpu::Atomic::Add(&kcount[N/2], 1);
                        if(power > tolerance)
                        {
                            Gpu::Atomic::Add(&ps_map[N/2], power);
                        }

                        break;
                    }

                    // If you've reached the largest bin but not been captured
                    else if(s > N/2)
                    { 
                        Print() << iv << "\n";
                        Print() << kmag << "\n";
                        Print() << kiso[s] << "," << kiso[s-1] << "\n";
                        Error("RandomField::print_power_spectrum Part of the spectrum isn't captured.");
                    }

                    // If you haven't found the right bin yet
                    else { continue; }
                }
            }
        });
    }

    ParallelAllReduce::Sum(kcount.data(), static_cast<int>(kcount.size()), ParallelContext::CommunicatorSub());
    ParallelAllReduce::Sum(ps_map.data(), static_cast<int>(ps_map.size()), ParallelContext::CommunicatorSub());

    // Print the power spectrum to a new file in data/
#pragma omp single
    for(int s=0; s<=N/2; s++)
    {
        power_spec_file.write_data_line({(kiso[s]+kiso[s+1])/2., (double)ps_map[s]/kcount[s]});
    }
}

// Finds statistical moment x of given MultiFab
inline Real RandomField::find_field_moment_x(MultiFab &field, const Vector<Real> mean, 
                                             const int moment, const int component)
{
    Real sum = 0.;
    const Real vol = std::pow(N, 3.);

    // Loop over the box to calculate moment x
    for (MFIter mfi(field); mfi.isValid(); ++mfi) 
    {
        Array4<Real> const& field_ptr = field.array(mfi);
        const Box& bx = mfi.fabbox();

        amrex::ParallelFor(bx, [=, &sum] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            sum += std::pow(field_ptr(i, j, k, component) - mean[component], moment);
        });
    }

    ParallelAllReduce::Sum(sum, ParallelContext::CommunicatorSub());
    //if (moment == 3) { Print() << "Components of skewness: ";
    //                   Print() << sum << ", " << sum/vol << "\n"; }

    // Normalise and return moment x
    if (sum == 0) { return 0; }
    else if(moment == 2) { return sqrt(sum/vol); }
    else { return sum/vol; }
}

// Calculates and prints requested moments (any between 1 and 4)
inline Vector<Real> RandomField::print_moment(MultiFab &field, const Vector<std::string> names,  
                                             const Vector<int> &moment_orders, SmallDataIO &file, 
                                             const int is_first_step)
{
    // Trap instance where the user requests too large a moment
    for(const auto moment : moment_orders)
    {
        if(moment > 4) 
        { 
            Error("RandomField::print_moment Chosen moment order has not been implemented");
        }
    }

    // Allocate arrays to store each moment
    const int nc = field.nComp();
    const Real vol = std::pow(N, 3.);
    Vector<Real> means(nc, 0.);
    Vector<Real> stdev(nc, 0.);
    Vector<Real> skew(nc, 0.);
    Vector<Real> kurt(nc, 0.);

    // Find iterators, which determine which moments are requested and their ordering
    Vector<int>::const_iterator start = moment_orders.begin();
    Vector<int>::const_iterator mean_itr = std::find(moment_orders.begin(), moment_orders.end(), 1);
    Vector<int>::const_iterator stdev_itr = std::find(moment_orders.begin(), moment_orders.end(), 2);
    Vector<int>::const_iterator skew_itr = std::find(moment_orders.begin(), moment_orders.end(), 3);
    Vector<int>::const_iterator kurt_itr = std::find(moment_orders.begin(), moment_orders.end(), 4);

    // Allocate vectors to store header line and data lines
    Vector<Real> data_to_print(nc * moment_orders.size(), 0.);
    Vector<std::string> headers(nc * moment_orders.size(), "");

    for (int comp = 0; comp < nc; comp++)
    {
        means[comp] = field.sum(comp)/vol;
        if(mean_itr != moment_orders.end())
        {
            assign_statistics_data(headers, names[comp]+"-mean", data_to_print, means, comp, nc,
                                    mean_itr, start, is_first_step);
        }

        if(moment_orders.back() != 1)
        {
            stdev[comp] = find_field_moment_x(field, means, 2, comp);
            if(stdev_itr != moment_orders.end())
            {
                assign_statistics_data(headers, names[comp]+"-stdev", data_to_print, stdev, comp, nc, 
                                        stdev_itr, start, is_first_step);
            }

            if(moment_orders.back() != 2)
            {
                skew[comp] = find_field_moment_x(field, means, 3, comp);
                skew[comp] /= std::pow(stdev[comp], 3.);

                if(skew_itr != moment_orders.end())
                {
                    assign_statistics_data(headers, names[comp]+"-skew", data_to_print, skew, comp, nc,
                                            skew_itr, start, is_first_step);
                }

                if(moment_orders.back() != 3)
                {
                    kurt[comp] = find_field_moment_x(field, means, 4, comp);
                    kurt[comp] /= std::pow(stdev[comp], 4.);

                    assign_statistics_data(headers, names[comp]+"-kurt", data_to_print, kurt, comp, nc,
                                            kurt_itr, start, is_first_step);
                }
            }
        }
    }

    // Write header and data lines
#pragma omp single
    if(is_first_step) 
    { 
        file.write_header_line(headers); 
    }

#pragma omp single
    file.write_time_data_line(data_to_print);
    
    return stdev;
}

inline void RandomField::derive(const MultiFab &source, MultiFab &out, int dcomp)
{
    BL_PROFILE("RandomField::derive");

    // Extract MultiFab ingredients from state
    BoxArray sba = source.boxArray();
    DistributionMapping sdm = source.DistributionMap();
    MultiFab hij_x(sba, sdm, 6, 0);
    hij_x.setVal(0.0);

    // Copy the spatial metric from the state
    Copy(hij_x, source, c_h11, lut[0][0], 1, 0);
    Copy(hij_x, source, c_h12, lut[0][1], 1, 0);
    Copy(hij_x, source, c_h13, lut[0][2], 1, 0);
    Copy(hij_x, source, c_h22, lut[1][1], 1, 0);
    Copy(hij_x, source, c_h23, lut[1][2], 1, 0);
    Copy(hij_x, source, c_h33, lut[2][2], 1, 0);

    // Undo the normalisation and BSSN-CPT conversion
    for (int l=0; l<3; l++) { hij_x.plus(-1., lut[l][l], 1); }
    hij_x.mult(1./norm);

    // Set up the problem domain in Fourier space
    // And impose that MPI ranks only slice along the i index (for Nyquist conditions)
    IntVect domain_low(0, 0, 0);
    IntVect k_domain_high(N/2, N-1, N-1);
    Box k_domain(domain_low, k_domain_high);
    Array< bool, AMREX_SPACEDIM > const &slicing{true, false, false};
    BoxArray kba = decompose(k_domain, ParallelContext::NProcsAll(), slicing);
    DistributionMapping kdm(kba);

    // Set up the arrays to store the Fourier data sets
    cMultiFab hs_k(kba, kdm, 2, 0);
    cMultiFab hij_k(kba, kdm, 6, 0);
    hs_k.setVal(0.0);
    hij_k.setVal(0.0);

    // Set up the FFT
    IntVect x_domain_high(N-1, N-1, N-1);
    Box x_domain(domain_low, x_domain_high);
    FFT::R2C<Real> tensor_fft(x_domain, FFT::Info().setBatchSize(hij_k.nComp()));

    // Perform the fft
    tensor_fft.forward(hij_x, hij_k);

    // Normalise the fft (fftw style)
    for(int comp = 0; comp < 6; comp++)
    {
        hij_k.mult(std::pow(N, -3.), comp, 1); 
    }

    // Loop to extract the Fourier-space mode functions
    for (MFIter mfi(hij_k); mfi.isValid(); ++mfi) 
    {
        const Box& bx = mfi.fabbox();

        // Make a pointer to the mode functions at this MF box
        Array4<GpuComplex<Real>> const& hs_ptr = hs_k.array(mfi);
        Array4<GpuComplex<Real>> const& hij_ptr = hij_k.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            IntVect iv{i, j, k};
            Vector<Real> mhat(3, 0.);
            Vector<Real> nhat(3, 0.);

            mhat = calculate_basis_vector(iv, 0);
            nhat = calculate_basis_vector(iv, 1);

            Real eplus = 0.;
            Real ecross = 0.;

            // Find basis tensors and do the Fourier trick
            for (int l=0; l<3; l++) for (int p=0; p<3; p++)
            {
                eplus = mhat[l]*mhat[p] - nhat[l]*nhat[p];
                ecross = mhat[l]*nhat[p] + nhat[l]*mhat[p];

                hs_ptr(i, j, k, 0) += (hij_ptr(i, j, k, lut[l][p]) * eplus)/std::sqrt(2.);
                hs_ptr(i, j, k, 1) += (hij_ptr(i, j, k, lut[l][p]) * ecross)/std::sqrt(2.);
            }
        });
    }

    apply_nyquist_conditions(hs_k);

    // Make a multifab to store config space mode functions
    // Need to use out to make these ingredients??
    BoxArray xba = out.boxArray();//(x_domain); //
    DistributionMapping xdm = out.DistributionMap();//(xba); //
    MultiFab hs_x(xba, xdm, 2, 0);
    hs_x.setVal(0.0);

    // Fourier transform
    FFT::R2C<Real> mode_function_fft(x_domain, FFT::Info().setBatchSize(hs_k.nComp()));
    mode_function_fft.backward(hs_k, hs_x);

    // Apply physical normalisation
    hs_x.mult(norm);

    for (MFIter mfi(hs_x); mfi.isValid(); ++mfi) 
    {
        Array4<Real> const& out_ptr = out.array(mfi);
        Array4<Real> const& hx_ptr = hs_x.array(mfi);

        const Box& bx = mfi.fabbox();
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const IntVect iv{i, j, k};
            bool in_ghost_index = is_ghost_index(iv);
            if(!in_ghost_index)
            {
                out_ptr(iv, dcomp) = hx_ptr(i, j, k, 0);
                out_ptr(iv, dcomp + 1) = hx_ptr(i, j, k, 1);
            }
        });
    }
}

// Main extraction routine
inline void RandomField::extract(const MultiFab &state, const std::string data_path, const Real dt,  
                                 const Real cur_time, const int restart_time, const int first_step)
{
    BL_PROFILE("RandomField::extract");

    // Extract MultiFab ingredients from state
    BoxArray sba = state.boxArray();
    DistributionMapping sdm = state.DistributionMap();
    MultiFab gij_x(sba, sdm, 6, 0);

    // 0: scalar field
    // 1: conformal factor
    MultiFab scalars_x(sba, sdm, 2, 0);

    // Copy the spatial metric from the state
    Copy(gij_x, state, c_h11, lut[0][0], 1, 0);
    Copy(gij_x, state, c_h12, lut[0][1], 1, 0);
    Copy(gij_x, state, c_h13, lut[0][2], 1, 0);
    Copy(gij_x, state, c_h22, lut[1][1], 1, 0);
    Copy(gij_x, state, c_h23, lut[1][2], 1, 0);
    Copy(gij_x, state, c_h33, lut[2][2], 1, 0);

    int m_c_phi = 0;
    int m_c_chi = 1;
    Copy(scalars_x, state, c_phi, m_c_phi, 1, 0);
    Copy(scalars_x, state, c_chi, m_c_chi, 1, 0);

    // Find background quantities needed to extract \cal R
    const int vol = std::pow(m_params.N_readin, 3);
    const double K_bar = state.sum(c_K)/vol;
    const double alpha_bar = state.sum(c_lapse)/vol;
    const double Pi_bar = state.sum(c_Pi)/vol;
    const double phi_bar = state.sum(c_phi)/vol;

    // Remove background from scalar field
    scalars_x.plus(-phi_bar, m_c_phi, 1);
    scalars_x.plus(-1., m_c_chi, 1);
    scalars_x.mult(1./norm);

    // Undo the normalisation and BSSN-CPT conversion
    for (int l=0; l<3; l++) { gij_x.plus(-1., lut[l][l], 1); }
    gij_x.mult(1./norm);

    // Set up the problem domain in Fourier space
    // And impose that MPI ranks only slice along the i index (for Nyquist conditions)
    IntVect domain_low(0, 0, 0);
    IntVect k_domain_high(N/2, N-1, N-1);
    Box k_domain(domain_low, k_domain_high);
    Array< bool, AMREX_SPACEDIM > const &slicing{true, false, false};
    BoxArray kba = decompose(k_domain, ParallelContext::NProcsAll(), slicing);
    DistributionMapping kdm(kba);

    // Set up the arrays to store the Fourier data sets
    cMultiFab hs_k(kba, kdm, 2, 0);
    cMultiFab gij_k(kba, kdm, 6, 0);
    cMultiFab scalars_k(kba, kdm, 2, 0);
    cMultiFab R_k(kba, kdm, 1, 0);

    hs_k.setVal(0.0);
    gij_k.setVal(0.0);
    scalars_k.setVal(0.0);
    R_k.setVal(0.0);

    // Set up the FFT
    IntVect x_domain_high(N-1, N-1, N-1);
    Box x_domain(domain_low, x_domain_high);
    FFT::R2C<Real> tensor_fft(x_domain, FFT::Info().setBatchSize(gij_k.nComp()));
    FFT::R2C<Real> scalar_fft(x_domain, FFT::Info().setBatchSize(scalars_k.nComp()));

    // Perform the fft
    tensor_fft.forward(gij_x, gij_k);
    scalar_fft.forward(scalars_x, scalars_k);

    // Normalise the fft (fftw style)
    for(int comp = 0; comp < 6; comp++) { gij_k.mult(std::pow(N, -3.), comp, 1); }
    for(int comp = 0; comp < 2; comp++) { scalars_k.mult(std::pow(N, -3.), comp, 1); }

    std::string Filename = "/cephfs/home/eaf49/GRTeclyn-workspace/Examples/ScalarField/Rk-check.dat";

    // Loop to extract the Fourier-space mode functions
    for (MFIter mfi(gij_k); mfi.isValid(); ++mfi) 
    {
        const Box& bx = mfi.fabbox();

        // Make a pointer to the mode functions at this MF box
        Array4<GpuComplex<Real>> const& hs_ptr = hs_k.array(mfi);
        Array4<GpuComplex<Real>> const& hij_ptr = gij_k.array(mfi);
        Array4<GpuComplex<Real>> const& scalars_ptr = scalars_k.array(mfi);
        Array4<GpuComplex<Real>> const& R_k_ptr = R_k.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            IntVect iv{i, j, k};
            
            Vector<Real> mhat(3, 0.);
            Vector<Real> nhat(3, 0.);

            mhat = calculate_basis_vector(iv, 0);
            nhat = calculate_basis_vector(iv, 1);

            Tensor<2, Real> eplus, ecross;

            // Find basis tensors and do the Fourier trick
            for (int l=0; l<3; l++) for (int p=0; p<3; p++)
            {
                eplus[l][p] = mhat[l]*mhat[p] - nhat[l]*nhat[p];
                ecross[l][p] = mhat[l]*nhat[p] + nhat[l]*mhat[p];

                hs_ptr(i, j, k, 0) += (hij_ptr(i, j, k, lut[l][p]) * eplus[l][p])/2.;
                hs_ptr(i, j, k, 1) += (hij_ptr(i, j, k, lut[l][p]) * ecross[l][p])/2.;
            }

            Tensor<2, GpuComplex<Real>> hij, hSV;
            GpuComplex<Real> hij_tr = 0.;
            GpuComplex<Real> hSV_tr = 0.;
            for (int l=0; l<3; l++) for (int p=0; p<3; p++)
            {
                hij[l][p] = hs_ptr(i, j, k, 0) * eplus[l][p] + hs_ptr(i, j, k, 1) * ecross[l][p];
                hSV[l][p] = hij_ptr(i, j, k, lut[l][p]) - hij[l][p];

                if(l==p) { hij_tr += hij[l][p]; hSV_tr += hSV[l][p]; }

                // if(i == 1 && j == 0 && k == 1)
                // {
                //     Print() << "In extract: \n";
                //     Print() << l << ", " << p << ": ";
                //     Print() << hij[l][p] << ", " << hSV[l][p] << "\n";
                // }
            }
            if(i == 1 && j == 2 && k == 1)
            {
                Print() << "Traces: ";
                Print() << hij_tr << ", " << hSV_tr << "\n";
            }

            if(m_params.scalar_init)
            {
                Vector<Real> iv_k(iv.begin(), iv.end());
                for(auto& k_comp : iv_k) { k_comp *= 2. * M_PI / m_params.L; }
                Real kmag = get_kmag(iv);
                GpuComplex<Real> Phi = 0;

                if(kmag == 0)
                {
                    R_k_ptr(i, j, k, 0) = GpuComplex<Real>{0., 0.};
                }

                else
                {
                    // converstion from chi, gamma_ij -> Phi
                    for(int l=0; l<3; l++) for(int p=0; p<3; p++)
                    {
                        Phi += (iv_k[l] * iv_k[p] * hSV[l][p])/std::pow(kmag, 2.);
                    }
                    Phi *= -1./48.;
                    Phi += 0.5 * (scalars_ptr(i, j, k, m_c_chi));

                    // calculate R_k
                    R_k_ptr(i, j, k, 0) = Phi - K_bar * scalars_ptr(i, j, k, m_c_phi) / alpha_bar / Pi_bar;
                }
            }
        });
    }

    apply_nyquist_conditions(hs_k);
    apply_nyquist_conditions(R_k);

    // Find the binned PS for each mode function and print to data/
    if((m_params.calc_binned_power_spectrum) 
	    && (static_cast<int>(cur_time/dt) % m_params.plot_int == 0))
    {
	    Print() << "Time step at print: " << static_cast<int>(std::round(cur_time/dt)) << "\n";
        std::string spec_path = make_subdirectory(data_path, "spectra", first_step);
        Vector<std::string> filenames(2, "");

        for(int comp = 0; comp < hs_k.nComp(); comp++)
        {
            filenames[comp] = spec_path+"spectrum-comp-"+std::to_string(comp)+"-time-";
            SmallDataIO spectrum_file(filenames[comp], dt, cur_time, restart_time, SmallDataIO::NEW, first_step, ".dat");
            print_power_spectrum(hs_k, spectrum_file, comp);
        }

        std::string filename = spec_path+"spectrum-Rk-time-";
        SmallDataIO spectrum_file(filename, dt, cur_time, restart_time, SmallDataIO::NEW, first_step, ".dat");
        print_power_spectrum(R_k, spectrum_file, 0);
    }

    // Find mode functions in configuration space if requested
    if(m_params.calc_higher_order_statistics)
    {
        // Make a multifab to store config space mode functions
        BoxArray xba(x_domain);
        DistributionMapping xdm(xba);
        MultiFab hs_x(xba, xdm, 2, 0);
        MultiFab R_x(xba, xdm, 1, 0);
        hs_x.setVal(0.0);
        R_x.setVal(0.0);

        // Fourier transform
        FFT::R2C<Real> tensor_mode_function_fft(x_domain, FFT::Info().setBatchSize(hs_k.nComp()));
        FFT::R2C<Real> scalar_mode_function_fft(x_domain, FFT::Info().setBatchSize(R_k.nComp()));
        tensor_mode_function_fft.backward(hs_k, hs_x);
        scalar_mode_function_fft.backward(R_k, R_x);

        // Apply physical normalisation
        hs_x.mult(norm);
        R_x.mult(norm);

        int output_comps = hs_x.nComp() + R_x.nComp();
        MultiFab out_MF(hs_x.boxArray(), hs_x.DistributionMap(), output_comps, 0);
        Copy(out_MF, R_x, 0, 0, R_k.nComp(), 0);
        Copy(out_MF, hs_x, 0, R_k.nComp(), hs_x.nComp(), 0);

        // Print mode functions if requested
        if(m_params.print_mode_functions == 1)
        {
            std::string mf_path = make_subdirectory(data_path, "mode-functions", first_step);
            std::string filename = mf_path+"mode-function-"+std::to_string(cur_time/dt);

            for (MFIter mfi(hs_x); mfi.isValid(); ++mfi) 
            {
                Array4<Real> const& hx_ptr = hs_x.array(mfi);
                const Box& bx = mfi.fabbox();

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    AllPrintToFile(filename) << i + N*(j + N*k) << ",";
                    for(int c=0; c<2; c++)
                    {
                        AllPrintToFile(filename).SetPrecision(14) << hx_ptr(i, j, k, c) << ",";
                    }
                    AllPrintToFile(filename) << "\n";
                });
            }
        }

        // Calculate and print field moments if requested
        Vector<Real> stdevs;
        if (m_params.calc_higher_order_statistics)
        {
            SmallDataIO stats_file(data_path+"field-statistics", dt, cur_time, 
                                    restart_time, SmallDataIO::APPEND, first_step, ".dat");

            if (!m_params.orders.empty())
            {
                Vector<std::string> names{"R","hplus","hcross"};
                stdevs = print_moment(out_MF, names, m_params.orders, stats_file, first_step);
            }
        }
        
        SmallDataIO ts_file(data_path+"tensor-scalar-ratio", dt, cur_time, 
                                    restart_time, SmallDataIO::APPEND, first_step, ".dat");
#pragma omp single
        if(first_step) 
    	{ 
            ts_file.write_header_line({"T/S ratio (plus)", "T/S ratio (cross)"}); 
        }
        
#pragma omp single
    	ts_file.write_time_data_line({stdevs[1] / stdevs[0], stdevs[2] / stdevs[0]}); 
    }
}

#endif /* RANDOMFIELD_IMPL_HPP_*/
