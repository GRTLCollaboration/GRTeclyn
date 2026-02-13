/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef WORMHOLEINITIALDATA_HPP_
#define WORMHOLEINITIALDATA_HPP_

#include "CCZ4StateVariables.hpp" // For c_chi, c_h11, etc.
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "VarsTools.hpp"
#include "simd.hpp"

//! Class which calculates the initial data for a Morris-Thorne wormhole
class WormholeInitialData
{
  public:
    //! Structure for the input parameters
    struct params_t
    {
        double throat_radius; // b_0
        double k_amplitude;   // Amplitude of the Gaussian perturbation
        double k_width;       // Width (sigma) of the Gaussian
        std::array<double, AMREX_SPACEDIM> center;
    };

    //! Constructor
    WormholeInitialData(params_t a_params, double a_dx)
        : m_params(a_params), m_dx(a_dx)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t>
    void compute(int i, int j, int k, amrex::Array4<data_t> cell) const
    {
        // 1. Get Coordinates
        // The grid coordinates (x,y,z) correspond to the proper distance l
        amrex::IntVect grid_index(i, j, k);
        Coordinates<data_t> coords(grid_index, m_dx, m_params.center);

        data_t x = coords.x;
        data_t y = coords.y;
        data_t z = coords.z;

        // Calculate l^2 (proper radial distance squared)
        // Add a tiny epsilon to avoid division by zero at the exact center
        data_t l2 = x * x + y * y + z * z + 1e-12;
        data_t l  = std::sqrt(l2);

        // 2. Geometric Parameters
        double b0    = m_params.throat_radius;
        double b0_sq = b0 * b0;

        // A = b_0^2 / l^2
        data_t A_term = b0_sq / l2;

        // 3. Conformal Factor
        // chi = (1 + b_0^2/l^2)^(-2/3)
        data_t chi = std::pow(1.0 + A_term, -2.0 / 3.0);

        // 4. Conformal Metric \tilde{gamma}_{ij}
        // Derived form: \tilde{gamma}_{ij} = (1+A)^{1/3} delta_{ij} -
        // A(1+A)^{-2/3} n_i n_j Where n_i = x_i / l

        data_t prefactor_delta = std::pow(1.0 + A_term, 1.0 / 3.0);
        data_t prefactor_n     = A_term * chi; // A * (1+A)^(-2/3)

        // Normal vector components n_i
        data_t nx = x / l;
        data_t ny = y / l;
        data_t nz = z / l;

        Tensor<2, data_t> h;

        // Diagonal
        h[0][0] = prefactor_delta - prefactor_n * nx * nx;
        h[1][1] = prefactor_delta - prefactor_n * ny * ny;
        h[2][2] = prefactor_delta - prefactor_n * nz * nz;

        // Off-diagonal
        h[0][1] = -prefactor_n * nx * ny;
        h[0][2] = -prefactor_n * nx * nz;
        h[1][2] = -prefactor_n * ny * nz;

        // Symmetric parts
        h[1][0] = h[0][1];
        h[2][0] = h[0][2];
        h[2][1] = h[1][2];

        // 5. Extrinsic Curvature
        // K = K_amp * exp(-l^2 / sigma^2)
        // A_ij (traceless part) = 0
        data_t K = m_params.k_amplitude *
                   exp(-l2 / (m_params.k_width * m_params.k_width));

        // 6. Lapse and Shift
        // alpha = 1, beta = 0
        data_t lapse = 1.0;

        // 7. Store variables
        cell(i, j, k, c_chi) = chi;

        // Metric components
        cell(i, j, k, c_h11) = h[0][0];
        cell(i, j, k, c_h12) = h[0][1];
        cell(i, j, k, c_h13) = h[0][2];
        cell(i, j, k, c_h22) = h[1][1];
        cell(i, j, k, c_h23) = h[1][2];
        cell(i, j, k, c_h33) = h[2][2];

        // Extrinsic curvature
        cell(i, j, k, c_K) = K;
        // Traceless A_ij is zero
        cell(i, j, k, c_A11) = 0.0;
        cell(i, j, k, c_A12) = 0.0;
        cell(i, j, k, c_A13) = 0.0;
        cell(i, j, k, c_A22) = 0.0;
        cell(i, j, k, c_A23) = 0.0;
        cell(i, j, k, c_A33) = 0.0;

        // Gauge variables
        cell(i, j, k, c_lapse)  = lapse;
        cell(i, j, k, c_shift1) = 0.0;
        cell(i, j, k, c_shift2) = 0.0;
        cell(i, j, k, c_shift3) = 0.0;
        cell(i, j, k, c_B1)     = 0.0;
        cell(i, j, k, c_B2)     = 0.0;
        cell(i, j, k, c_B3)     = 0.0;

        // Theta (CCZ4 damping)
        cell(i, j, k, c_Theta) = 0.0;
    }

  protected:
    WormholeInitialData::params_t m_params;
    double m_dx;
};

#endif /* WORMHOLEINITIALDATA_HPP_ */