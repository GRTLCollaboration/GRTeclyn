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

inline bool RandomField::is_ghost_index(IntVect vector)
{
    bool ret = false;
    for(int d=0; d<3; d++) 
    { 
        if(vector[d] < 0 || vector[d] > N-1) { ret = true; }
    }
    return ret;
}

inline std::string RandomField::make_subdirectory(std::string base, std::string dir, int is_first_step)
{
    std::string new_path = base+dir+"/";
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

inline void RandomField::assign_statistics_data(Vector<std::string> &header_storage, const std::string name, 
                            Vector<Real> &data_storage, const Array1D<Real, 0, 1> data, 
                            const int component, const auto itr, const auto start, const int is_first_step)
{
    int loc = component + 2*(itr - start);
    if(is_first_step) 
    { 
        header_storage[loc] =  name+std::to_string(component); 
    }
    data_storage[loc] = data(component);
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
        value *= 0.5 * (1.0 - tanh(Dt * (kmag - ks))); 
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

inline bool RandomField::is_independent_draw(IntVect iv)
{
    bool val = true;
    int i = iv[0];
    int j = iv[1];
    int k = iv[2];

    if((i==0 || i==N/2) && (j==0 || j==N/2) && (k==0 || k== N/2)) { val = false; }
    else if (i==0 || i==N/2) 
    {
        if((k>N/2 && j==N/2) || (k==0 && j>N/2) || (k>N/2 && j==0) || (k==N/2 && j>N/2))
        {
            val = false;
        }
        else if(j > N/2)
        {
            val = false;
        }
    }

    return val;
}

inline int RandomField::find_processor(const IntVect iv, const BoxArray &ba, Vector<int> pm)
{
    int proc;
    bool found = false;

    for(int p=0; p<pm.size(); p++)
    {
        Box this_box = ba[p];
        if(iv >= this_box.smallEnd() && iv < this_box.bigEnd())
        {
            proc = p;
            found = true;
        }
    }

    if(found == false)
    {
        Print() << iv << "\n";
        Print() << ba;
        Error("RandomField::find_processor Processor not found.");
    }

    return proc;
}

inline void RandomField::apply_nyquist_conditions(cMultiFab &field, const BoxArray &ba, const DistributionMapping &dm)
{
    //REDESIGN so that the loop is over only one MPI rank/the MPI memory is shared

    Vector<int> pm = dm.ProcessorMap();
    const int nc = field.nComp();

    for (MFIter mfi(field); mfi.isValid(); ++mfi) 
    {
        Array4<GpuComplex<Real>> const& field_ptr = field.array(mfi);
        const Box& bx = mfi.fabbox();

        /*amrex::ParallelFor(bx, [=, &mfi, &field] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if ((i==0 || i==N/2) && (j==0 || j==N/2) && (k==0 || k== N/2))
            {
                for(int comp = 0; comp < nc; comp++)
                {
                    GpuComplex<Real> temp(field(i, j, k, comp).real(), 0.);
                    field(i, j, k, comp) = temp;
                }
            }
            else 
            {
                IntVect iv{i, j, k};
                IntVect iv_invt{i, invert_index(j), invert_index(k)};
                IntVect iv_flp{i, invert_index(j), flip_index(k)};

                int invt_proc = find_processor(iv_invt, ba, pm);
                int flip_proc = find_processor(iv_flp, ba, pm);
                Print(0) << mfi.index() << ", " << invt_proc << ", " << flip_proc << "\n";
            }
        });*/
    }

    //Error();



    /*Print() << "BA: " << ba << "\n";
    
    Print() << "DM size: " << pm.size() << "\n";
    for(int p=0; p < pm.size(); p++)
    {
        Print() << p << ", " << pm[p] << "\n";
        Print() << ba[pm[p]] << "\n";
    }*/

    //Print() << plane_data[0][0] << "\n";

    /*Array4<GpuComplex<Real>> values;

    for (MFIter mfi(field); mfi.isValid(); ++mfi) 
    {
        Array4<GpuComplex<Real>> const& field_ptr = field.array(mfi);

        for(int i=0; i<=N/2; i++) for(int j=0; j<N; j++) for(int k=0; k<N; k++) for(int comp = 0; comp < field.nComp(); comp++)
        {
            if (i==0 || i==N/2) 
            {
                if((k>N/2 && j==N/2) || (k==0 && j>N/2) || (k>N/2 && j==0) || (k==N/2 && j>N/2))
                {
                    GpuComplex<Real> temp(field_ptr(i, invert_index(j), invert_index(k), comp).real(), 
                                            -field_ptr(i, invert_index(j), invert_index(k), comp).imag());
                    values(i, j, k, comp) = temp;
                }
                else if(j > N/2)
                {
                    GpuComplex<Real> temp(field_ptr(i, invert_index(j), flip_index(k), comp).real(), 
                                            -field_ptr(i, invert_index(j), flip_index(k), comp).imag());
                    values(i, j, k, comp) = temp;
                }
                else { values(i, j, k, comp) = 0.; }
            }
        }
    }*/

    // Bcast(field(i,j,k,comp))

    // Nyquist node condition
    /*if ((i==0 || i==N/2) && (j==0 || j==N/2) && (k==0 || k== N/2))
    {
        for(int comp = 0; comp < field.nComp(); comp++)
        {
            GpuComplex<Real> temp(field(i, j, k, comp).real(), 0.);
            field(i, j, k, comp) = temp;
        }
    }*/

    // Nyquist axis condition
    /*if (i==0 || i==N/2) 
    {
        if((k>N/2 && j==N/2) || (k==0 && j>N/2) || (k>N/2 && j==0) || (k==N/2 && j>N/2))
        {
            for(int comp = 0; comp < field.nComp(); comp++) 
            {
                GpuComplex<Real> temp(field(i, invert_index(j), invert_index(k), comp).real(), 
                                        //-field(i, invert_index(j), invert_index(k), comp).imag());
                field(i, j, k, comp) = temp;
            }
        }
        /*else if(j > N/2)
        {
            for(int comp = 0; comp < field.nComp(); comp++) 
            {
                GpuComplex<Real> temp(field(i, invert_index(j), flip_index(k), comp).real(), 
                                        -field(i, invert_index(j), flip_index(k), comp).imag());
                field(i, j, k, comp) = temp;
            }
        }
    }*/
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

    std::string Filename = "/nfs/st01/hpc-gr-epss/eaf49/GRTeclyn-dump/hs-k-init";

    Print() << "MFs initialised, starting initial data loop.\n";

    GpuComplex<Real> zero(0., 0.);

    const IntVect start1{0, 0, 0};
    const IntVect end1{0, N-1, N-1};
    Box nyq_plane_1(start1, end1);

    const IntVect start2{0, 0, 0};
    const IntVect end2{0, N-1, N-1};
    Box nyq_plane_2(start1, end1);

    BaseFab<GpuComplex<Real>> nyq_bf_1(nyq_plane_1, 6);
    BaseFab<GpuComplex<Real>> nyq_bf_2(nyq_plane_1, 6);
    
    
    /*Dim3 const& start{0, 0, 0};
    Dim3 const& end{N-1, N-1, N-1};
    Array4<GpuComplex<Real>> nyquist_plane_values(&zero, start, end, 6);*/

    //Vector<GpuComplex<Real>> nyquist_plane_point(6, zero);
    //Vector<Vector<GpuComplex<Real>>> nyquist_plane_values(2.*std::pow(N-1, 2.), nyquist_plane_point);

    /*const IntVect start1{0, 0, 0};
    const IntVect end1{0, N-1, N-1};
    Box plane_dim1(start1, end1);
    BoxArray nyquist_array1(plane_dim1);
    DistributionMapping nyquist_dm(nyquist_array1, 1);
    cMultiFab nyquist_plane_1(nyquist_array1, nyquist_dm, 6, 0);*/
    //Vector<int> pm = nyquist_dm.ProcessorMap();

    //Print() << nyquist_array1;
    /*Print() << nyquist_dm.ProcessorMap().size() << "\n";
    for(int proc : nyquist_dm.ProcessorMap())
    {
        Print() << proc << "\n";
    }*/
    //Print() << nyquist_dm.ProcessorMap();
    /*Error();*/

    //std::map<std::tuple<int, int, int>, std::tuple<int, GpuComplex<Real>>> nyquist_plane_values;

    //const auto &plane_arrays = nyquist_plane_1.arrays();

    // Loop to create Fourier-space tensor object
    //MFIter::allowMultipleMFIters(true);
    for (MFIter mfi(hs_k); mfi.isValid(); ++mfi) 
    {
        // Make a pointer to the mode functions at this MF box
        Array4<GpuComplex<Real>> const& hs_ptr = hs_k.array(mfi);
        Array4<GpuComplex<Real>> const& hij_ptr = hij_k.array(mfi);

        Array4<GpuComplex<Real>> const& As_ptr = As_k.array(mfi);
        Array4<GpuComplex<Real>> const& Aij_ptr = Aij_k.array(mfi);

        Array4<GpuComplex<Real>> const& nyq_array_1 = nyq_bf_1.array();
        Array4<GpuComplex<Real>> const& nyq_array_2 = nyq_bf_2.array();

        //Array4<GpuComplex<Real>> const& nyq_ptr = nyquist_plane_1.array(mfi);

        //MFIter mfin(nyquist_plane_1);
        //Array4<GpuComplex<Real>> const& nyquist_ptr_1 = nyquist_plane_1.array(mfin);

        //const auto &ptr = nyquist_plane_1.arrays();
        //Array4<GpuComplex<Real>> const& nyquist_ptr_1 = ptr[1];
        //Print() << nyquist_ptr_1.size();// << ", " << nyquist_ptr_1.nComp() << "\n";
        //std::cout << ptr[0] << "\n";

        const Box& bx = mfi.fabbox();
        Print() << "Tensor min index: " << bx.smallEnd() << "\n";
        Print() << "Tensor max index: " << bx.bigEnd() << "\n";

        // Loop to create mode functions then hij(k)
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            IntVect iv = {i, j, k};
            bool independent_draw = is_independent_draw(iv);

            if(independent_draw)
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

                if (i==0) 
                {
                    for(int comp=0; comp < 6; comp++)
                    {
                        nyq_array_1(i, j, k, comp) = hij_ptr(i, j, k, comp);
                    }
                    //CellData<GpuComplex<Real>> nyq_cell = plane_arrays[0].cellData(i, j, k);
                    //std::cout << nyq_cell.nComp() << "\n";
                    //for(int comp=0; comp < 6; comp++)
                    //{
                        //nyq_cell[comp] = 0.;
                        //nyquist_ptr(i, j, k, comp) = hij_ptr(i, j, k, comp);
                    //}
                    /*for (MFIter mfin(hs_k); mfin.isValid(); ++mfin) 
                    {
                        Array4<GpuComplex<Real>> const& nyquist_ptr_1 = nyquist_plane_1.array(mfin);
                        for(int comp=0; comp < 6; comp++)
                        {
                            nyquist_ptr_1(i, j, k, comp) = hij_ptr(i, j, k, comp);
                        }*/

                        
                        /*auto iterator = ps_map.find(std::tuple<i, j, k>);
                        auto& [component, field] = iterator->second;
                        Gpu::Atomic::Exch(&component, comp);
                        Gpu::Atomic::Exch(&field, hij_ptr(i, j, k, comp));

                        //nyquist_plane_values[k + N * j][comp] = hij_ptr(i, j, k, comp);

                        //std::cout << iv << ": " << comp << ", " << nyquist_plane_values(i, j, k, comp) << "\n";// = zero;//hij_ptr(i, j, k, comp); 
                    }*/
                }
                else if (i == N/2)
                {
                    for(int comp=0; comp < 6; comp++)
                    {
                        nyq_array_2(i, j, k, comp) = hij_ptr(i, j, k, comp);
                    }
                }
            }

            /*IntVect iv{i, j, k};
            bool in_ghost_index = is_ghost_index(iv);
            if(!in_ghost_index)
            {
                for(int s=0; s<2; s++) 
                { 
                    PrintToFile(Filename, 0) << iv << ": ";
                    PrintToFile(Filename, 0) << hs_ptr(i, j, k, s).real() << ","; 
                    PrintToFile(Filename, 0) << hs_ptr(i, j, k, s).imag() << ",";
                }
                PrintToFile(Filename, 0) << "\n";
            }*/
        });
    }

    //ParallelDescriptor::Bcast(nyquist_plane_values.dataPtr(), nyquist_plane_values.size());
    //ParallelDescriptor::Bcast(nyquist_plane_values[0].data(), nyquist_plane_values[0].size());

    //apply_nyquist_conditions(hs_k);
    apply_nyquist_conditions(hij_k, sba, sdm);
    //Error();
    //apply_nyquist_conditions(Aij_k);

    Print() << "Starting Fourier transform.\n";

    for(int fcomp = 0; fcomp < hij_k.nComp(); fcomp++)
    {
        cMultiFab hij_k_slice(hij_k, make_alias, fcomp, 1);
        MultiFab hij_x_slice(hij_x, make_alias, fcomp, 1);
        random_field_fft.backward(hij_k_slice, hij_x_slice);

        cMultiFab Aij_k_slice(Aij_k, make_alias, fcomp, 1);
        MultiFab Aij_x_slice(Aij_x, make_alias, fcomp, 1);
        random_field_fft.backward(Aij_k_slice, Aij_x_slice);
    }

    Print() << "Starting normalisation.\n";

    hij_x.mult(norm);
    Aij_x.mult(norm);

    for (int l=0; l<3; l++) { hij_x.plus(1., lut[l][l], 1); }
    Aij_x.mult(-0.5);

    Print() << "Putting ICs on the grid.\n";

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

    Print() << "Initialisation ended.\n";
}

/****
    Extraction routines
****/

inline Real RandomField::find_field_moment_x(MultiFab &field, Array1D<Real, 0, 1> mean, 
                                                int moment, int component)
{
    Real sum = 0.;
    const Real vol = std::pow(N, 3.);

    for (MFIter mfi(field); mfi.isValid(); ++mfi) 
    {
        Array4<Real> const& field_ptr = field.array(mfi);
        const Box& bx = mfi.fabbox();

        amrex::ParallelFor(bx, [=, &sum] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            sum += std::pow(field_ptr(i, j, k, component) - mean(component), moment);
        });
    }
    if(moment == 2) { return sqrt(sum/vol); }
    else { return sum/vol; }
}

inline void RandomField::print_tensor_moment(MultiFab &field, const Vector<int> &moment_orders, 
                                    SmallDataIO &file, const int is_first_step)
{
    for(const auto moment : moment_orders)
    {
        if(moment > 4) 
        { 
            Error("RandomField::print_tensor_moment Chosen moment order has not been implemented");
        }
    }

    const Real vol = std::pow(N, 3.);
    Array1D<Real, 0, 1> means = {0., 0.};
    Array1D<Real, 0, 1> stdev = {0., 0.};
    Array1D<Real, 0, 1> skew = {0., 0.};
    Array1D<Real, 0, 1> kurt = {0., 0.};

    const auto start = moment_orders.begin();
    const auto mean_itr = std::find(moment_orders.begin(), moment_orders.end(), 1);
    const auto stdev_itr = std::find(moment_orders.begin(), moment_orders.end(), 2);
    const auto skew_itr = std::find(moment_orders.begin(), moment_orders.end(), 3);
    const auto kurt_itr = std::find(moment_orders.begin(), moment_orders.end(), 4);

    Vector<Real> data_to_print(2 * moment_orders.size(), 0.);
    Vector<std::string> headers(2 * moment_orders.size(), "");

    for(int comp = 0; comp < 2; comp++)
    {
        means(comp) = field.sum(comp)/vol;
        if(mean_itr != moment_orders.end())
        {
            assign_statistics_data(headers, "Mean", data_to_print, means, comp, 
                                    mean_itr, start, is_first_step);
        }

        if(moment_orders.back() != 1)
        {
            stdev(comp) = find_field_moment_x(field, means, 2, comp);
            if(stdev_itr != moment_orders.end())
            {
                assign_statistics_data(headers, "Stdev", data_to_print, stdev, comp, 
                                        stdev_itr, start, is_first_step);
            }

            if(moment_orders.back() != 2)
            {
                skew(comp) = find_field_moment_x(field, means, 3, comp);
                skew(comp) /= std::pow(stdev(comp), 3.);

                if(skew_itr != moment_orders.end())
                {
                    assign_statistics_data(headers, "Skew", data_to_print, skew, comp,
                                            skew_itr, start, is_first_step);
                }

                if(moment_orders.back() != 3)
                {
                    kurt(comp) = find_field_moment_x(field, means, 4, comp);
                    kurt(comp) /= std::pow(stdev(comp), 4.);

                    assign_statistics_data(headers, "Kurt", data_to_print, kurt, comp,
                                            kurt_itr, start, is_first_step);
                }
            }
        }
    }

    if(is_first_step) { file.write_header_line(headers); }
    file.write_time_data_line(data_to_print);
}

inline void RandomField::print_power_spectrum(cMultiFab &field_array, SmallDataIO &power_spec_file, int component)
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

    // Set up isotropic k axis and PS map
    double dk_to_bin = (double)m_params.bin_number/((double)N/2);
    double kmag = 0.;
    Vector<Real> kiso(N/2+1, 0.);
    std::map<int, std::tuple<Real, Real, int>> ps_map;

    for (int s=0; s<=N/2; s++) { kiso[s] = s*dkiso; ps_map[s] = std::make_tuple(kiso[s], 0., 0); }

    // Loop to bin the power spectrum at each point
    MFIter::allowMultipleMFIters(true); // Needed to pass the map to the ParallelFor loop
    for (MFIter mfi(field_array); mfi.isValid(); ++mfi) 
    {
        // Make a pointer to the mode function at this MF box
        Array4<GpuComplex<Real>> const& field_ptr = field_array.array(mfi);
        const Box& bx = mfi.fabbox();

        amrex::ParallelFor(bx, [=, &ps_map] AMREX_GPU_DEVICE (int i, int J, int K) noexcept
        {
            IntVect iv{i, J, K};

            // Check to see if you're in a ghost cell
            bool in_ghost_index = is_ghost_index(iv);
            if(!in_ghost_index)
            {
                int j = invert_index(J);
                int k = invert_index(K);
                double kmag = std::sqrt(i*i + j*j + k*k) * 2 * M_PI / m_params.L;

                // make sure you're still in the domain
                if(kmag > kiso_max) 
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
                        auto iterator = ps_map.find(s-1);
                        auto& [iso, sum, kcount] = iterator->second;

                        Real power = (std::pow(field_ptr(i, j, k, component).real(), 2.0) 
                                                    + std::pow(field_ptr(i, j, k, component).imag(), 2.0));
                        Gpu::Atomic::Add(&sum, power);
                        Gpu::Atomic::Add(&kcount, 1);

                        break;
                    }

                    // If you're at the largest bin
                    else if(s == N/2)
                    { 
                        auto iterator = ps_map.find(N/2);
                        auto& [iso, sum, kcount] = iterator->second;

                        Real power = (std::pow(field_ptr(i, j, k, component).real(), 2.0) 
                                                    + std::pow(field_ptr(i, j, k, component).imag(), 2.0));
                        Gpu::Atomic::Add(&sum, power);
                        Gpu::Atomic::Add(&kcount, 1);

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

    // Print the power spectrum to a new file in data/
    for(int s=0; s<=N/2; s++)
    {
        auto [isotropic_k, power, count] = ps_map[s];
        power_spec_file.write_data_line({isotropic_k, (double)power/count});
    }
}

inline void RandomField::extract(MultiFab &state, std::string data_path, Real dt, Real cur_time, int restart_time, int first_step)
{
    BL_PROFILE("RandomField::init");

    // Extract MultiFab ingredients from state
    BoxArray sba = state.boxArray();
    DistributionMapping sdm = state.DistributionMap();
    MultiFab hij_x(sba, sdm, 6, 0);

    // Copy the spatial metric from the state
    Copy(hij_x, state, c_h11, lut[0][0], 1, 0);
    Copy(hij_x, state, c_h12, lut[0][1], 1, 0);
    Copy(hij_x, state, c_h13, lut[0][2], 1, 0);
    Copy(hij_x, state, c_h22, lut[1][1], 1, 0);
    Copy(hij_x, state, c_h23, lut[1][2], 1, 0);
    Copy(hij_x, state, c_h33, lut[2][2], 1, 0);

    // Undo the normalisation and delta function 
    for (int l=0; l<3; l++) { hij_x.plus(-1., lut[l][l], 1); }
    hij_x.mult(1./norm);

    // Set up the problem domain and MF ingredients (Fourier space)
    IntVect domain_low(0, 0, 0);
    IntVect domain_high(N-1, N-1, N-1);
    Box domain(domain_low, domain_high);
    FFT::R2C<Real> random_field_fft(domain);
    auto const& [kba, kdm] = random_field_fft.getSpectralDataLayout();

    // Set up the arrays to store the Fourier data sets
    cMultiFab hs_k(kba, kdm, 2, 0);
    cMultiFab hij_k(kba, kdm, 6, 0);

    // Perform the fft
    for(int fcomp = 0; fcomp < hij_x.nComp(); fcomp++)
    {
        cMultiFab hij_k_slice(hij_k, make_alias, fcomp, 1);
        MultiFab hij_x_slice(hij_x, make_alias, fcomp, 1);
        random_field_fft.forward(hij_x_slice, hij_k_slice);
    }

    // Normalise the fft (fftw style)
    for(int comp = 0; comp < 6; comp++)
    {
        hij_k.mult(std::pow(N, -3.), comp, 1); 
    }

    std::string Filename = "/nfs/st01/hpc-gr-epss/eaf49/GRTeclyn-dump/hs-k-extr";

    // Loop to extract the Fourier-space mode functions
    for (MFIter mfi(hij_k); mfi.isValid(); ++mfi) 
    {
        // Make a pointer to the mode functions at this MF box
        Array4<GpuComplex<Real>> const& hs_ptr = hs_k.array(mfi);
        Array4<GpuComplex<Real>> const& hij_ptr = hij_k.array(mfi);
        const Box& bx = mfi.fabbox();

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            IntVect iv{i, j, k};
            Vector<Real> mhat(3, 0.);
            Vector<Real> nhat(3, 0.);

            mhat = calculate_basis_vector(i, j, k, 0);
            nhat = calculate_basis_vector(i, j, k, 1);

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

            /*bool in_ghost_index = is_ghost_index(iv);
            if(!in_ghost_index)
            {
                for(int s=0; s<2; s++) 
                { 
                    //PrintToFile(Filename, 0) << iv << ": ";
                    PrintToFile(Filename, 0).SetPrecision(12) << hs_ptr(i, j, k, s).real() << ","; 
                    PrintToFile(Filename, 0).SetPrecision(12) << hs_ptr(i, j, k, s).imag() << ",";
                }
                PrintToFile(Filename, 0) << "\n";
            }*/
        });

        // Apply Nyquist line and plane conditions
        /*amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            apply_nyquist_conditions(i, j, k, hs_ptr);
        });*/
    }

    //apply_nyquist_conditions(hs_k);

    // Find the binned PS for each mode function and print to data/
    if(m_params.calc_binned_power_spectrum) 
    {
        for(int comp = 0; comp < hs_k.nComp(); comp++)
        {
            std::string spec_path = make_subdirectory(data_path, "spectra", first_step);

            SmallDataIO spectrum_file(spec_path+"spectrum-comp-"+std::to_string(comp)+"-time-", 
                                        dt, cur_time, restart_time, SmallDataIO::NEW, first_step, ".dat");
            if(first_step) { spectrum_file.write_header_line({"k", "power"}, ""); }
            print_power_spectrum(hs_k, spectrum_file, comp);
        }
    }

    if(m_params.calc_config_space_mode_fns || m_params.calc_higher_order_statistics)
    {
        // Make a multifab to store config space mode functions
        BoxArray xba(domain);
        DistributionMapping xdm(xba);
        MultiFab hs_x(xba, xdm, 2, 0);

        for(int fcomp = 0; fcomp < hs_x.nComp(); fcomp++)
        {
            cMultiFab hs_k_slice(hs_k, make_alias, fcomp, 1);
            MultiFab hs_x_slice(hs_x, make_alias, fcomp, 1);
            random_field_fft.backward(hs_k_slice, hs_x_slice);
        }

        hs_x.mult(norm);

        if(m_params.calc_config_space_mode_fns)
        {
            std::string mf_path = make_subdirectory(data_path, "mode-functions", first_step);

            SmallDataIO mode_function_file(mf_path+"mode-function-", dt, cur_time, 
                                            restart_time, SmallDataIO::NEW, first_step, ".dat");
            
            if(first_step) { mode_function_file.write_header_line({"plus field", "cross field"}, ""); }
            for (MFIter mfi(hs_x); mfi.isValid(); ++mfi) 
            {
                Array4<Real> const& hx_ptr = hs_x.array(mfi);
                const Box& bx = mfi.fabbox();

                amrex::ParallelFor(bx, [=, &mode_function_file] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    mode_function_file.write_data_line({hx_ptr(i, j, k, 0), hx_ptr(i, j, k, 1)});
                });
            }
        }

        if (m_params.calc_higher_order_statistics)
        {
            SmallDataIO stats_file(data_path+"field-statistics", dt, cur_time, 
                                    restart_time, SmallDataIO::APPEND, first_step, ".dat");

            if (!m_params.orders.empty())
            {
                print_tensor_moment(hs_x, m_params.orders, stats_file, first_step);
            }
        }
    }
}

#endif /* RANDOMFIELD_IMPL_HPP_*/
