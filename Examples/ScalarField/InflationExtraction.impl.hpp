/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef INFLATIONEXTRACTION_IMPL_HPP_
#define INFLATIONEXTRACTION_IMPL_HPP_

/* Helper functions */

// Makes subdirectories in data/
inline std::string InflationExtraction::make_subdirectory(const std::string base, 
                                                          const std::string dir, 
                                                          const int is_m_first_step) const
{
    std::string new_path = base+"../"+dir+"/";
    if(is_m_first_step)
    {
        if (FilesystemTools::directory_exists(base)) 
        { 
            FilesystemTools::mkdir_recursive(new_path); 
        }
        else 
        { 
            Print() << "Directory creation failed for " << new_path << "\n";
            Error("InflationExtraction::extract Data directory has not been created."); 
        }
    }
    return new_path;
}

// Creates a custom data file layout 
inline void InflationExtraction::assign_statistics_data(Vector<std::string> &header_storage, const std::string name, 
                            Vector<Real> &data_storage, const Vector<Real> data, const int component, const int num_comps,
                            const Vector<int>::const_iterator itr, const Vector<int>::const_iterator start, 
                            const int is_m_first_step)
{
    int loc = component + num_comps*(itr - start);
    if(is_m_first_step) 
    { 
        header_storage[loc] =  name; 
    }
    data_storage[loc] = data[component];
}

// Calculates and prints the power spectrum
inline void InflationExtraction::print_power_spectrum(const cMultiFab &field_array, 
                                                      SmallDataIO &power_spec_file, 
                                                      const int component)
{ 
    // Set up the isotropic k axis bounds
    double kiso_max = std::sqrt(3.) * m_params.N * M_PI / m_params.L;
    double dkiso = sqrt(3.) * 2. * M_PI / m_params.L;

    // check the stepping along the diagonal is consistent
    if (kiso_max/dkiso - m_params.N/2 > InflationUtils::tolerance)
    {
        Error("RandomField::print_power_spectrum "
              "Isotropic k axis is too large.");
    }
    // check you aren't sampling above the max sampling rate
    else if (m_params.bin_number > kiso_max/dkiso)
    {
        Error("RandomField::print_power_spectrum "
              "Bin number is too large.");
    }
    // check your bin number isn't greater than the max resolvable bins
    else if(m_params.bin_number > m_params.N/2)
    {
        Error("RandomField::print_power_spectrum "
              "bin number must be less than N/2.");
    }

    // Set up isotropic k axis and PS map
    double dk_to_bin = m_params.bin_number / (m_params.N/2.);
    double kmag = 0.;
    Vector<Real> kiso(m_params.N / 2 + 1, 0.);

    Vector<Real> ps_map(m_params.bin_number + 1, 0.);
    Vector<int> kcount(m_params.bin_number + 1, 0);

    for (int s=0; s<=m_params.N/2; s++) { kiso[s] = s*dkiso; }

    // Needed to pass the map to the ParallelFor loop
    MFIter::allowMultipleMFIters(true); 

    // Loop to bin the power spectrum at each point
    const auto& field_arrs = field_array.arrays();
    ParallelFor(field_array, [=, &ps_map, &kcount]
                AMREX_GPU_DEVICE (int bx, int i, int j, int k)
        {
            // Check to see if you're in a ghost cell
            IntVect iv{i, j, k};
            Real kmag = m_params.get_kmag(iv);

            // make sure you're still in the domain
            if(kmag - kiso_max > InflationUtils::tolerance) 
            { 
                Print() << iv << "\n";
                Print() << kmag << "," << kiso_max << "\n";
                Error("RandomField::print_power_spectrum "
                      "Found magnitude larger than (N/2,N/2,N/2)."); 
            }

            // Loop over the isotropic axis
            for (int s=1; s<=m_params.N/2; s++) 
            {
                // If smaller than the smallest bin
                if(kmag < kiso[0])
                {
                    Print() << iv << "\n";
                    Error("RandomField::print_power_spectrum "
                          "kmag below the kiso domain.");
                }

                // If you're larger than the largest bin
                else if(kmag - kiso[m_params.N/2] > InflationUtils::tolerance)
                {
                    Print() << iv << "\n";
                    Error("RandomField::print_power_spectrum "
                          "kmag above the kiso domain.");
                }

                // If you're somewhere in the middle
                else if ((kmag < kiso[s] && kmag >= kiso[(s-1)]) 
                        || kmag == kiso[m_params.N/2]) 
                {
                    Real power = (std::pow(field_arrs[bx](i, j, k, component).real(), 2.0) 
                                + std::pow(field_arrs[bx](i, j, k, component).imag(), 2.0));
                    
                    int comp = (kmag == kiso[m_params.N/2]) ? m_params.N/2 : s - 1;
                    Gpu::Atomic::Add(&kcount[comp], 1);
                    if(power > InflationUtils::tolerance)
                    {
                        Gpu::Atomic::Add(&ps_map[comp], power);   
                    }

                    break;
                }

                // If you've reached the largest bin but not been captured
                else if(s > m_params.N/2)
                { 
                    Print() << iv << "\n";
                    Print() << kmag << "\n";
                    Print() << kiso[s] << "," << kiso[s-1] << "\n";
                    Error("RandomField::print_power_spectrum "
                          "Part of the spectrum isn't captured.");
                }

                // If you haven't found the right bin yet
                else { continue; }
            }
        });

    ParallelAllReduce::Sum(kcount.data(), static_cast<int>(kcount.size()), 
                                          ParallelContext::CommunicatorSub());
    ParallelAllReduce::Sum(ps_map.data(), static_cast<int>(ps_map.size()), 
                                          ParallelContext::CommunicatorSub());

    // Print the power spectrum to a new file in data/
#pragma omp single
    for(int s = 0; s <= m_params.N/2; s++)
    {
        power_spec_file.write_data_line({(kiso[s]+kiso[s+1])/2., ps_map[s]/kcount[s]});
    }
}

// Finds statistical moment x of given MultiFab
inline Real InflationExtraction::calculate_field_moment_x(const MultiFab &field, const Vector<Real> mean, 
                                                          const int moment, const int component)
{
    Real sum = 0.;
    const Real vol = std::pow(m_params.N, 3.);

    const auto& field_arrs = field.arrays();

    ParallelFor(field, [=, &sum] AMREX_GPU_DEVICE
                (int bx, int i, int j, int k)
        {
            sum += std::pow(field_arrs[bx](i, j, k, component) - mean[component], moment);
        });

    ParallelAllReduce::Sum(sum, ParallelContext::CommunicatorSub());

    // Normalise and return moment x
    if (sum == 0) { return 0; }
    else if(moment == 2) { return sqrt(sum/vol); }
    else { return sum/vol; }
}

// Calculates and prints requested moments (any between 1 and 4)
inline Vector<Real> InflationExtraction::print_moment(MultiFab &field, const Vector<std::string> names,  
                                             const Vector<int> &moment_orders, SmallDataIO &file, 
                                             const int is_first_step)
{
    // Trap instance where the user requests too large a moment
    for(const auto moment : moment_orders)
    {
        if(moment > 4) 
        { 
            Error("InflationExtraction::print_moment "
                  "Chosen moment order has not been implemented");
        }
    }

    // Allocate arrays to store each moment
    const int nc = field.nComp();
    const Real vol = std::pow(m_params.N, 3.);
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
            stdev[comp] = calculate_field_moment_x(field, means, 2, comp);
            if(stdev_itr != moment_orders.end())
            {
                assign_statistics_data(headers, names[comp]+"-stdev", data_to_print, stdev, comp, nc, 
                                        stdev_itr, start, is_first_step);
            }

            if(moment_orders.back() != 2)
            {
                skew[comp] = calculate_field_moment_x(field, means, 3, comp);
                skew[comp] /= std::pow(stdev[comp], 3.);

                if(skew_itr != moment_orders.end())
                {
                    assign_statistics_data(headers, names[comp]+"-skew", data_to_print, skew, comp, nc,
                                            skew_itr, start, is_first_step);
                }

                if(moment_orders.back() != 3)
                {
                    kurt[comp] = calculate_field_moment_x(field, means, 4, comp);
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

/* Main functions */

// Extract R and hs in configuration space from the BSSN variables
inline void InflationExtraction::extract_hs_and_R(MultiFab &hs, MultiFab &R, 
                                                  const MultiFab &state, const bool print_spec = false)
{
    // Extract MultiFab ingredients from state
    BoxArray sba = state.boxArray();
    DistributionMapping sdm = state.DistributionMap();
    if (sba != hs.boxArray() 
        || sdm != hs.DistributionMap())
    {
        Error("InflationExtraction::extract_hs_and_R "
              "source and output BA or SDM do not match");
    }

    // 0: scalar field
    // 1: conformal factor
    MultiFab scalars_x(sba, sdm, 2, 0);
    MultiFab gij_x(sba, sdm, 6, 0);

    // Copy the spatial metric from the state
    Copy(gij_x, state, c_h11, InflationUtils::lut[0][0], 1, 0);
    Copy(gij_x, state, c_h12, InflationUtils::lut[0][1], 1, 0);
    Copy(gij_x, state, c_h13, InflationUtils::lut[0][2], 1, 0);
    Copy(gij_x, state, c_h22, InflationUtils::lut[1][1], 1, 0);
    Copy(gij_x, state, c_h23, InflationUtils::lut[1][2], 1, 0);
    Copy(gij_x, state, c_h33, InflationUtils::lut[2][2], 1, 0);

    int m_c_phi = 0;
    int m_c_chi = 1;
    Copy(scalars_x, state, c_phi, m_c_phi, 1, 0);
    Copy(scalars_x, state, c_chi, m_c_chi, 1, 0);

    // Find background quantities needed to extract \cal R
    const int vol = std::pow(m_params.N, 3);
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
    for (int l=0; l<3; l++) { gij_x.plus(-1., InflationUtils::lut[l][l], 1); }
    gij_x.mult(1./norm);

    // Set up the problem domain in Fourier space
    // And impose that MPI ranks only slice along the i index (for Nyquist conditions)
    IntVect domain_low(0, 0, 0);
    IntVect k_domain_high(m_params.N/2, m_params.N-1, m_params.N-1);
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
    IntVect x_domain_high(m_params.N-1, m_params.N-1, m_params.N-1);
    Box x_domain(domain_low, x_domain_high);
    FFT::R2C<Real> tensor_fft(x_domain, FFT::Info().setBatchSize(gij_k.nComp()));
    FFT::R2C<Real> scalar_fft(x_domain, FFT::Info().setBatchSize(scalars_k.nComp()));

    // Perform the fft
    tensor_fft.forward(gij_x, gij_k);
    scalar_fft.forward(scalars_x, scalars_k);

    // Normalise the fft (fftw style)
    for(int comp = 0; comp < 6; comp++) { gij_k.mult(std::pow(m_params.N, -3.), comp, 1); }
    for(int comp = 0; comp < 2; comp++) { scalars_k.mult(std::pow(m_params.N, -3.), comp, 1); }

    // Set variables to store the maximum trace 
    // of the scalar and tensor components
    Real hij_tr_max = 0.;
    Real hSV_tr_max = 0.;

    const auto& hs_arrs = hs_k.arrays();
    const auto& gij_arrs = gij_k.arrays();
    const auto& scalars_arrs = scalars_k.arrays();
    const auto& R_k_arrs = R_k.arrays();

    ParallelFor(gij_k, [=, &hij_tr_max, &hSV_tr_max]
                AMREX_GPU_DEVICE (int bx, int i, int j, int k)
        {
            IntVect iv{i, j, k};
            if (iv != IntVect{0, 0, 0})
            {
                Tensor<2, Real> eplus = m_params.calculate_polarisation_tensor(iv, 0);
                Tensor<2, Real> ecross = m_params.calculate_polarisation_tensor(iv, 1);

                // Find basis tensors and do the Fourier trick
                for (int l=0; l<3; l++) for (int p=0; p<3; p++)
                {
                    hs_arrs[bx](i, j, k, 0) += (gij_arrs[bx](i, j, k, InflationUtils::lut[l][p])
                                                 * eplus[l][p])/2.;
                    hs_arrs[bx](i, j, k, 1) += (gij_arrs[bx](i, j, k, InflationUtils::lut[l][p])
                                                 * ecross[l][p])/2.;
                }

                // Calculate the TT and scalar-(vector) components of the 
                // metric, by reconstructing hij and subtracting it from \tilde{gamma}_ij
                Tensor<2, GpuComplex<Real>> hij, hSV;
                GpuComplex<Real> hij_tr = 0.;
                GpuComplex<Real> hSV_tr = 0.;
                for (int l=0; l<3; l++) for (int p=0; p<3; p++)
                {
                    hij[l][p] = (hs_arrs[bx](i, j, k, 0) * eplus[l][p] 
                                + hs_arrs[bx](i, j, k, 1) * ecross[l][p]);
                    hSV[l][p] = gij_arrs[bx](i, j, k, InflationUtils::lut[l][p]) - hij[l][p];

                    // Find the traces of these components, as a diagnostic
                    if(l==p) { hij_tr += hij[l][p]; hSV_tr += hSV[l][p]; }
                }

                hij_tr_max = max(sqrt(pow(hij_tr.real(), 2.) + pow(hij_tr.imag(), 2.)),
                                 hij_tr_max);
                hSV_tr_max = max(sqrt(pow(hSV_tr.real(), 2.) + pow(hSV_tr.imag(), 2.)),
                                 hSV_tr_max);

                // Confirm hij is trace-free in Fourier space
                if (std::abs(hij_tr_max) > InflationUtils::tolerance)
                {
                    Print() << iv << ": " << hij_tr_max << "\n";
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
                    Real kmag = m_params.get_kmag(iv);
                    GpuComplex<Real> Phi = 0;

                    // converstion from chi and gamma_ij -> Phi
                    for(int l=0; l<3; l++) for(int p=0; p<3; p++)
                    {
                        Phi += (iv_k[l] * iv_k[p] * hSV[l][p])/std::pow(kmag, 2.);
                    }
                    Phi *= 1./4.;
                    Phi += 0.5 * (scalars_arrs[bx](i, j, k, m_c_chi));

                    // Combine the above to find R(k)
                    R_k_arrs[bx](i, j, k, 0) = Phi - (K_bar * scalars_arrs[bx](i, j, k, m_c_phi) 
                                                        / alpha_bar / Pi_bar);
                }
            }
        });

    // Output the max traces of the tensor components as a diagnostic
    SmallDataIO trace_file(m_data_path+"tensor-traces", m_dt, m_cur_time, 
                                m_restart_time, SmallDataIO::APPEND, m_first_step, ".dat");
    if(m_first_step) 
    { 
        trace_file.write_header_line({"hij trace max", "hSV trace max"}); 
    }
    trace_file.write_time_data_line({hij_tr_max, hSV_tr_max}); 

    // Prepare to IFT the polarisation fields and R field
    m_params.apply_nyquist_conditions(hs_k);
    m_params.apply_nyquist_conditions(R_k);

    // Find the binned PS for each mode function and print to data/
    if ((print_spec) && (static_cast<int>(m_cur_time/m_dt) 
                         % m_params.plot_int == 0))
    {
	    Print() << "Time step at print: " << static_cast<int>(std::round(m_cur_time/m_dt)) << "\n";
        std::string spec_path = make_subdirectory(m_data_path, "spectra", m_first_step);
        Vector<std::string> filenames(2, "");

        for(int comp = 0; comp < hs_k.nComp(); comp++)
        {
            filenames[comp] = spec_path+"spectrum-comp-"+std::to_string(comp)+"-time-";
            SmallDataIO spectrum_file(filenames[comp], m_dt, m_cur_time, m_restart_time, SmallDataIO::NEW, m_first_step, ".dat");
            print_power_spectrum(hs_k, spectrum_file, comp);
        }

        std::string filename = spec_path+"spectrum-Rk-time-";
        SmallDataIO spectrum_file(filename, m_dt, m_cur_time, m_restart_time, SmallDataIO::NEW, m_first_step, ".dat");
        print_power_spectrum(R_k, spectrum_file, 0);
    }

    // Fourier transform
    tensor_fft.backward(hs_k, hs);
    scalar_fft.backward(R_k, R);

    // Apply physical normalisation
    hs.mult(norm);
    R.mult(norm);
}

// Put R and hs into plotfiles
inline void InflationExtraction::derive(const MultiFab &source, MultiFab &out, const int dcomp)
{
    BL_PROFILE("InflationExtraction::derive");

    // Make a multifab to store config space mode functions
    BoxArray oba = out.boxArray();
    DistributionMapping odm = out.DistributionMap();
    MultiFab hs_x(oba, odm, 2, 0);
    MultiFab R_x(oba, odm, 1, 0);
    hs_x.setVal(0.0);
    R_x.setVal(0.0);

    extract_hs_and_R(hs_x, R_x, source);

    const auto& hs_arrs = hs_x.arrays();
    const auto& R_arrs = R_x.arrays();
    const auto& out_arrs = out.arrays();

    ParallelFor(hs_x, [=] AMREX_GPU_DEVICE
                (int bx, int i, int j, int k)
        {
            const IntVect iv{i, j, k};
            out_arrs[bx](iv, dcomp) = R_arrs[bx](i, j, k);
            out_arrs[bx](iv, dcomp + 1) = hs_arrs[bx](i, j, k, 0);
            out_arrs[bx](iv, dcomp + 2) = hs_arrs[bx](i, j, k, 1);
        });
}

// Find spectrum and higher-order statistics on R and hs
inline void InflationExtraction::extract(const MultiFab &state)
{
    BL_PROFILE("InflationExtraction::extract");

    // Make a multifab to store config space mode functions
    BoxArray sba = state.boxArray();
    DistributionMapping sdm = state.DistributionMap();
    MultiFab hs_x(sba, sdm, 2, 0);
    MultiFab R_x(sba, sdm, 1, 0);
    hs_x.setVal(0.0);
    R_x.setVal(0.0);

    extract_hs_and_R(hs_x, R_x, state, true);

    // Find mode functions in configuration space if requested,
    // and find the statistics (orders 1-4) of the polarisation fields 
    // and the R field. 
    // And print the tensor-to-scalar ratio if requested.
//     int output_comps = hs_x.nComp() + R_x.nComp();
//     MultiFab out_MF(hs_x.boxArray(), hs_x.DistributionMap(), output_comps, 0);
//     Copy(out_MF, R_x, 0, 0, R_x.nComp(), 0);
//     Copy(out_MF, hs_x, 0, R_x.nComp(), hs_x.nComp(), 0);

//     // Print mode functions if requested
//     if(m_params.print_mode_functions == 1)
//     {
//         std::string mf_path = make_subdirectory(m_data_path, "mode-functions", m_first_step);
//         std::string filename = mf_path+"mode-function-"+std::to_string(m_cur_time/m_dt);

//         for (MFIter mfi(hs_x); mfi.isValid(); ++mfi) 
//         {
//             Array4<Real> const& hx_ptr = hs_x.array(mfi);
//             const Box& bx = mfi.fabbox();

//             amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
//             {
//                 AllPrintToFile(filename) << i + m_params.N*(j + m_params.N*k) << ",";
//                 for(int c=0; c<2; c++)
//                 {
//                     AllPrintToFile(filename).SetPrecision(14) << hx_ptr(i, j, k, c) << ",";
//                 }
//                 AllPrintToFile(filename) << "\n";
//             });
//         }
//     }

//     // Calculate and print field moments
//     Vector<Real> stdevs;
//     SmallDataIO stats_file(m_data_path+"field-statistics", m_dt, m_cur_time, 
//                             m_restart_time, SmallDataIO::APPEND, m_first_step, ".dat");

//     if (!m_params.orders.empty())
//     {
//         Vector<std::string> names{"R","hplus","hcross"};
//         stdevs = print_moment(out_MF, names, m_params.orders, stats_file, m_first_step);
//     }
    
//     // Calculate and print tensor to scalar ratio (integrated PS)
//     SmallDataIO ts_file(m_data_path+"tensor-scalar-ratio", m_dt, m_cur_time, 
//                                 m_restart_time, SmallDataIO::APPEND, m_first_step, ".dat");
// #pragma omp single
//     if(m_first_step) 
//     { 
//         ts_file.write_header_line({"T/S ratio (plus)", "T/S ratio (cross)"}); 
//     }
    
// #pragma omp single
    // ts_file.write_time_data_line({stdevs[1] / stdevs[0], stdevs[2] / stdevs[0]}); 
}

#endif /* INFLATIONEXTRACTION_IMPL_HPP_ */