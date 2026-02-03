/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(WEYL4_HPP_)
#error "This file should only be included through Weyl4.hpp"
#endif

#ifndef WEYL4_IMPL_HPP_
#define WEYL4_IMPL_HPP_

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
Weyl4::operator()(int ix, int iy, int iz,
                  const amrex::Array4<amrex::Real> &weyl_scalars,
                  const amrex::Array4<amrex::Real const> &state) const
{
    // We do not want to amend the cell data values, so use const CCZ4Vars
    const amrex::CellData<const amrex::Real> &state_cell_data =
        state.cellData(ix, iy, iz);
    CCZ4Vars vars(state_cell_data);

    // we need d1 chi, K, h, A... this just gets all of them
    const CCZ4D1Vars d1(ix, iy, iz, state, m_deriv);
    // we only need d2 of chi and h
    const Tensor<2, amrex::Real> d2_chi =
        m_deriv.diff2(ix, iy, iz, state, c_chi);
    const Tensor<4, amrex::Real> d2_h =
        m_deriv.diff2_tensor(ix, iy, iz, state, c_h11);
    const auto h_UU  = CCZ4Geometry::compute_inverse_metric(vars);
    const auto chris = CCZ4Geometry::compute_christoffel(d1, h_UU);

    // Get the coordinates
    const Coordinates coords(amrex::IntVect(ix, iy, iz), m_dx, m_center);

    // Compute the spatial volume element
    const auto epsilon3_LUU = compute_epsilon3_LUU(vars, h_UU);

    // Compute the E and B fields
    EBFields_t ebfields =
        compute_EB_fields(vars, d1, d2_chi, d2_h, epsilon3_LUU, h_UU, chris);

    // work out the Newman Penrose scalar
    weyl_scalar_t out = compute_Weyl4(ebfields, vars, h_UU, coords);

    // store the result
    weyl_scalars(ix, iy, iz, m_out_comp)     = out.Real;
    weyl_scalars(ix, iy, iz, m_out_comp + 1) = out.Im;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor<3, amrex::Real>
Weyl4::compute_epsilon3_LUU(const CCZ4Vars &vars,
                            const Tensor<2, amrex::Real> &h_UU) const
{
    // raised normal vector, NB index 3 is time
    Tensor<1, amrex::Real, 4> n_U;
    n_U[3] = 1. / vars.lapse();
    FOR (i)
    {
        n_U[i] = -vars.shift(i) / vars.lapse();
    }

    // 4D levi civita symbol and 3D levi civita tensor in LLL and LUU form
    const auto epsilon4 = TensorAlgebra::epsilon4D();
    Tensor<3, amrex::Real> epsilon3_LLL;
    Tensor<3, amrex::Real> epsilon3_LUU;

    // Projection of antisymmentric Tensor onto hypersurface - see 8.3.17,
    // Alcubierre
    FOR (i, j, k)
    {
        epsilon3_LLL[i][j][k] = 0.0;
        epsilon3_LUU[i][j][k] = 0.0;
    }
    // projection of 4-antisymetric tensor to 3-tensor on hypersurface
    // note last index contracted as per footnote 86 pg 290 Alcubierre
    FOR (i, j, k)
    {
        for (int l = 0; l < 4; ++l)
        {
            epsilon3_LLL[i][j][k] += n_U[l] * epsilon4[i][j][k][l] *
                                     vars.lapse() /
                                     (vars.chi() * std::sqrt(vars.chi()));
        }
    }
    // rasing indices
    FOR (i, j, k)
    {
        FOR (m, n)
        {
            epsilon3_LUU[i][j][k] += epsilon3_LLL[i][m][n] * h_UU[m][j] *
                                     vars.chi() * h_UU[n][k] * vars.chi();
        }
    }

    return epsilon3_LUU;
}

// Calculation of E and B fields, using tetrads from gr-qc/0104063
// BSSN expressions from Alcubierre book
// CCZ4 expressions calculated by MR and checked with TF see:
// https://www.overleaf.com/read/tvqjbyhvqqtp
AMREX_GPU_DEVICE AMREX_FORCE_INLINE EBFields_t Weyl4::compute_EB_fields(
    const CCZ4Vars &vars, const CCZ4D1Vars &d1,
    const Tensor<2, amrex::Real> &d2_chi, const Tensor<4, amrex::Real> &d2_h,
    const Tensor<3, amrex::Real> &epsilon3_LUU,
    const Tensor<2, amrex::Real> &h_UU, const chris_t &chris) const
{
    EBFields_t out;

    // Extrinsic curvature
    Tensor<2, amrex::Real> K_tensor;
    Tensor<3, amrex::Real> d1_K_tensor;
    Tensor<3, amrex::Real> covariant_deriv_K_tensor;

    // Compute inverse, Christoffel symbols, Ricci tensor and Z terms
    // Note that unlike in CCZ4 equations we want R_ij + 0.5(D_iZ_j + D_jZ_i)
    // rather than R_ij + D_iZ_j + D_jZ_i hence use compute_ricci_Z_general
    double dZ_coeff        = (m_formulation == CCZ4RHS<>::USE_CCZ4) ? 1. : 0.;
    auto ricci_and_Z_terms = CCZ4Geometry::compute_ricci_Z_general(
        vars, d1, d2_chi, d2_h, h_UU, chris, dZ_coeff);

    // Compute full spatial Christoffel symbols
    const Tensor<3, amrex::Real> chris_phys =
        CCZ4Geometry::compute_phys_chris(d1.chi(), vars, h_UU, chris.ULL);

    // Extrinsic curvature and corresponding covariant and partial derivatives
    FOR (i, j)
    {
        K_tensor[i][j] = vars.A(i, j) / vars.chi() +
                         1. / 3. * (vars.h(i, j) * vars.K()) / vars.chi();

        FOR (k)
        {
            d1_K_tensor[i][j][k] =
                d1.A(i, j, k) / vars.chi() -
                d1.chi(k) / vars.chi() * K_tensor[i][j] +
                1. / 3. * d1.h(i, j, k) * vars.K() / vars.chi() +
                1. / 3. * vars.h(i, j) * d1.K(k) / vars.chi();
        }
    }
    // covariant derivative of K
    FOR (i, j, k)
    {
        covariant_deriv_K_tensor[i][j][k] = d1_K_tensor[i][j][k];

        FOR (l)
        {
            covariant_deriv_K_tensor[i][j][k] +=
                -chris_phys[l][k][i] * K_tensor[l][j] -
                chris_phys[l][k][j] * K_tensor[i][l];
        }
    }

    // Use 'K-Theta' in CCZ4. Just 'K' in BSSN. Not a mistake, this is not to
    // confuse with the typical 'K-2*Theta' that appears in the CCZ4 equations
    amrex::Real K_minus_Theta = vars.K();
    if (m_formulation == CCZ4RHS<>::USE_CCZ4)
    {
        K_minus_Theta -= vars.Theta();
    }

    // Calculate electric and magnetic fields
    FOR (i, j)
    {
        out.E[i][j] = 0.0;
        out.B[i][j] = 0.0;
    }

    FOR (i, j)
    {
        out.E[i][j] +=
            ricci_and_Z_terms.LL[i][j] + K_minus_Theta * K_tensor[i][j];

        FOR (k, l)
        {
            out.E[i][j] +=
                -K_tensor[i][k] * K_tensor[l][j] * h_UU[k][l] * vars.chi();

            out.B[i][j] +=
                epsilon3_LUU[i][k][l] * covariant_deriv_K_tensor[l][j][k];
        }
    }

    // For CCZ4, explicit symmetrization appears naturally;
    // For BSSN, only extra matter terms appear in the original expression, but
    // assuming the Momentum constraints are satisfied, we can make the
    // expression explicitly symmetric, which we enforce below
    // (see Alcubierre chapter 8.3, from eq. 8.3.15 onwards)
    TensorAlgebra::make_symmetric(out.B);

    // For CCZ4, Eij is explicitly trace-free;
    // For BSSN, only extra matter terms appear in the original expression, but
    // assuming the Hamiltonian constraint is satisfied, we can make the
    // expression explicitly trace free, which we enforce below
    // (see Alcubierre chapter 8.3, from eq. 8.3.15 onwards)
    CCZ4Geometry::make_trace_free(out.E, vars, h_UU);

    return out;
}

// Calculation of the Weyl4 scalar
AMREX_GPU_DEVICE AMREX_FORCE_INLINE weyl_scalar_t Weyl4::compute_Weyl4(
    const EBFields_t &ebfields, const CCZ4Vars &vars,
    const Tensor<2, amrex::Real> &h_UU, const Coordinates &coords) const
{
    weyl_scalar_t out;

    // Calculate the tetrads
    const Tetrad_t tetrad = compute_null_tetrad(vars, h_UU, coords);

    // Projection of Electric and magnetic field components using tetrads
    out.Real = 0.0;
    out.Im   = 0.0;
    FOR (i, j)
    {
        out.Real += 0.5 * (ebfields.E[i][j] * (tetrad.w[i] * tetrad.w[j] -
                                               tetrad.v[i] * tetrad.v[j]) -
                           2.0 * ebfields.B[i][j] * tetrad.w[i] * tetrad.v[j]);
        out.Im   += 0.5 * (ebfields.B[i][j] * (-tetrad.w[i] * tetrad.w[j] +
                                             tetrad.v[i] * tetrad.v[j]) -
                         2.0 * ebfields.E[i][j] * tetrad.w[i] * tetrad.v[j]);
    }

    return out;
}

// Calculation of the null tetrad
// Defintions from gr-qc/0104063
// "The Lazarus project: A pragmatic approach to binary black hole evolutions",
// Baker et al.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tetrad_t Weyl4::compute_null_tetrad(
    const CCZ4Vars &vars, const Tensor<2, amrex::Real> &h_UU,
    const Coordinates &coords) const
{
    Tetrad_t out;

    // compute coords
    const amrex::Real x = coords.x;
    const amrex::Real y = coords.y;
    const amrex::Real z = coords.z;

    // the alternating levi civita symbol
    const Tensor<3, double> epsilon = TensorAlgebra::epsilon();

    // calculate the tetrad
    out.u[0] = x;
    out.u[1] = y;
    out.u[2] = z;

    out.v[0] = -y;
    out.v[1] = x;
    out.v[2] = 0.0;

    out.w[0] = 0.0;
    out.w[1] = 0.0;
    out.w[2] = 0.0;

    FOR (i, j, k, m)
    {
        out.w[i] += 1. / std::sqrt(vars.chi()) * h_UU[i][j] * epsilon[j][k][m] *
                    out.v[k] * out.u[m];
    }

    // Gram Schmitt orthonormalisation
    // Choice of orthonormalisaion to avoid frame-dragging
    amrex::Real omega_11 = 0.0;
    FOR (i, j)
    {
        omega_11 += out.v[i] * out.v[j] * vars.h(i, j) / vars.chi();
    }
    FOR (i)
    {
        out.v[i] = out.v[i] / std::sqrt(omega_11);
    }

    amrex::Real omega_12 = 0.0;
    FOR (i, j)
    {
        omega_12 += out.v[i] * out.u[j] * vars.h(i, j) / vars.chi();
    }
    FOR (i)
    {
        out.u[i] += -omega_12 * out.v[i];
    }

    amrex::Real omega_22 = 0.0;
    FOR (i, j)
    {
        omega_22 += out.u[i] * out.u[j] * vars.h(i, j) / vars.chi();
    }
    FOR (i)
    {
        out.u[i] = out.u[i] / std::sqrt(omega_22);
    }

    amrex::Real omega_13 = 0.0;
    amrex::Real omega_23 = 0.0;
    FOR (i, j)
    {
        omega_13 += out.v[i] * out.w[j] * vars.h(i, j) / vars.chi();
        omega_23 += out.u[i] * out.w[j] * vars.h(i, j) / vars.chi();
    }
    FOR (i)
    {
        out.w[i] += -(omega_13 * out.v[i] + omega_23 * out.u[i]);
    }

    amrex::Real omega_33 = 0.0;
    FOR (i, j)
    {
        omega_33 += out.w[i] * out.w[j] * vars.h(i, j) / vars.chi();
    }
    FOR (i)
    {
        out.w[i] = out.w[i] / std::sqrt(omega_33);
    }

    return out;
}

AMREX_FORCE_INLINE
void Weyl4::set_up(int a_state_index)
{
    int num_ghosts = 2; // no advection terms so only need 2 ghost cells

    auto &derive_lst     = amrex::AmrLevel::get_derive_lst();
    const auto &desc_lst = amrex::AmrLevel::get_desc_lst();

    // Add Weyl4 to the derive list
    derive_lst.add(
        name, amrex::IndexType::TheCellType(),
        static_cast<int>(var_names.size()), var_names, Weyl4::compute_mf,
        [=](const amrex::Box &box) { return amrex::grow(box, 2); },
        &amrex::cell_quartic_interp);

    // We need all of the CCZ4 variables to calculate Weyl4
    // (except B but easier to keep it in)
    derive_lst.addComponent(name, desc_lst, a_state_index, 0, NUM_CCZ4_VARS);
}

AMREX_FORCE_INLINE
void Weyl4::compute_mf(amrex::MultiFab &out_mf, int out_comp, int ncomp,
                       const amrex::MultiFab &src_mf,
                       const amrex::Geometry &geomdata, amrex::Real /*time*/,
                       const int * /*bcrec*/, int /*level*/)
{
    const auto &out_arrays = out_mf.arrays();
    const auto &src_arrays = src_mf.const_arrays();

    GRParmParse pp;
    std::array<double, AMREX_SPACEDIM> center{};
    int formulation = 0;
    pp.get("extraction_center", center);
    pp.get("formulation", formulation);

    Weyl4 weyl4(center, geomdata.CellSize(0), out_comp, formulation);
    amrex::ParallelFor(
        out_mf, out_mf.nGrowVect(),
        [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz) noexcept
        { weyl4(ix, iy, iz, out_arrays[box_no], src_arrays[box_no]); });
}

#endif /* WEYL4_HPP_ */
