// Makes subdirectories in data/
inline std::string InflationExtraction::make_subdirectory(const std::string base, const std::string dir, const int is_first_step)
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
inline void InflationExtraction::assign_statistics_data(Vector<std::string> &header_storage, const std::string name, 
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

// Main extraction routine
inline void InflationExtraction::extract(const MultiFab &state, const std::string data_path, const Real dt,  
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
    const double chi_bar = state.sum(c_chi)/vol;

    // Remove background from scalar field
    scalars_x.plus(-phi_bar, m_c_phi, 1);
    scalars_x.plus(-chi_bar, m_c_chi, 1);
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

    // Set variables to store the maximum trace 
    // of the scalar and tensor components
    Real hij_tr_max = 0.;
    Real hSV_tr_max = 0.;

    // Loop to extract the Fourier-space mode functions
    for (MFIter mfi(gij_k); mfi.isValid(); ++mfi) 
    {
        const Box& bx = mfi.fabbox();
        Array4<GpuComplex<Real>> const& hs_ptr = hs_k.array(mfi);
        Array4<GpuComplex<Real>> const& hij_ptr = gij_k.array(mfi);
        Array4<GpuComplex<Real>> const& scalars_ptr = scalars_k.array(mfi);
        Array4<GpuComplex<Real>> const& R_k_ptr = R_k.array(mfi);

        amrex::ParallelFor(bx, [=, &hij_tr_max, &hSV_tr_max] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            IntVect iv{i, j, k};
            
            Vector<Real> mhat = calculate_basis_vector(iv, 0);
            Vector<Real> nhat = calculate_basis_vector(iv, 1);
            Tensor<2, Real> eplus, ecross;

            // Find basis tensors and do the Fourier trick
            for (int l=0; l<3; l++) for (int p=0; p<3; p++)
            {
                eplus[l][p] = mhat[l]*mhat[p] - nhat[l]*nhat[p];
                ecross[l][p] = mhat[l]*nhat[p] + nhat[l]*mhat[p];

                hs_ptr(i, j, k, 0) += (hij_ptr(i, j, k, lut[l][p]) * eplus[l][p])/2.;
                hs_ptr(i, j, k, 1) += (hij_ptr(i, j, k, lut[l][p]) * ecross[l][p])/2.;
            }

            if (m_params.alpha != 0) { Test_polarisation_tensor_orthonorm(iv, eplus, ecross); }

            // Calculate the TT and scalar-(vector) components of the 
            // metric, by reconstructing hij and subtracting it from \tilde{gamma}_ij
            Tensor<2, GpuComplex<Real>> hij, hSV;
            GpuComplex<Real> hij_tr = 0.;
            GpuComplex<Real> hSV_tr = 0.;
            for (int l=0; l<3; l++) for (int p=0; p<3; p++)
            {
                hij[l][p] = hs_ptr(i, j, k, 0) * eplus[l][p] + hs_ptr(i, j, k, 1) * ecross[l][p];
                hSV[l][p] = hij_ptr(i, j, k, lut[l][p]) - hij[l][p];

                // Find the traces of these components, as a diagnostic
                if(l==p) { hij_tr += hij[l][p]; hSV_tr += hSV[l][p]; }
            }
            
            Real hij_tr_mag = sqrt(pow(hij_tr.real(), 2.) + pow(hij_tr.imag(), 2.));
            Real hSV_tr_mag = sqrt(pow(hSV_tr.real(), 2.) + pow(hSV_tr.imag(), 2.));

            if (hij_tr_mag > hij_tr_max) { hij_tr_max = hij_tr_mag; }
            if (hSV_tr_mag > hSV_tr_max) { hSV_tr_max = hSV_tr_mag; }

            // Confirm hij is trace-free in Fourier space
            if (std::abs(hij_tr_mag) > tolerance)
            {
                Print() << iv << ": " << hij_tr_mag << "\n";
                Error("hij trace magnitude too large in extraction");
            }

            // Extract R according to the scheme detailed in 
            // Appendix B (Eq. B1) of arxiv:2502.06783, using hSV as the 
            // spatial metric instead of \tilde{gamma}_ij
            if(m_params.scalar_init)
            {
                // Find the unitful k vector
                Vector<Real> iv_k(iv.begin(), iv.end());
                for(auto& k_comp : iv_k) { k_comp *= 2. * M_PI / m_params.L; }
                Real kmag = get_kmag(iv);
                GpuComplex<Real> Phi = 0;

                // Set the zero mode
                if(kmag == 0)
                {
                    R_k_ptr(i, j, k, 0) = GpuComplex<Real>{0., 0.};
                }

                else
                {
                    // converstion from chi and gamma_ij -> Phi
                    for(int l=0; l<3; l++) for(int p=0; p<3; p++)
                    {
                        Phi += (iv_k[l] * iv_k[p] * hSV[l][p])/std::pow(kmag, 2.);
                    }
                    Phi *= 1./4.;
                    Phi += 0.5 * (scalars_ptr(i, j, k, m_c_chi));

                    // Combine the above to find R(k)
                    R_k_ptr(i, j, k, 0) = Phi - K_bar * scalars_ptr(i, j, k, m_c_phi) / alpha_bar / Pi_bar;
                }
            }
        });
    }

    // Output the max traces of the tensor components as a diagnostic
    SmallDataIO trace_file(data_path+"tensor-traces", dt, cur_time, 
                                restart_time, SmallDataIO::APPEND, first_step, ".dat");
    if(first_step) 
    { 
        trace_file.write_header_line({"hij trace max", "hSV trace max"}); 
    }
    trace_file.write_time_data_line({hij_tr_max, hSV_tr_max}); 

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

    // Find mode functions in configuration space if requested,
    // and find the statistics (orders 1-4) of the polarisation fields 
    // and the R field. 
    // And print the tensor-to-scalar ratio if requested.
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

        /* Print() << "Max tensor polarisations: " << hs_x.max(0) << ", " << hs_x.max(1) << "\n";
        Print() << "R max and min bounds: " << R_x.max(0) << ", " << R_x.min(0) << "\n"; */

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