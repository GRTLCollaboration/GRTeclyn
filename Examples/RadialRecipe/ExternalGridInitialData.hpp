/* ExternalGridInitialData
 * Reads a .gridinit binary file produced by the GRTresna bridge
 * and trilinear-interpolates it onto the AMReX grid.
 *
 * Supports both V1 (isotropic dx) and V2 (per-axis dx) formats:
 *   V1: "dx <scalar>"       → dx_x = dx_y = dx_z = scalar
 *   V2: "dx <dx_x> <dx_y> <dx_z>"
 *
 * Binary body: float64[nz][ny][nx][ncomp]   (C-order)
 *
 * The heavy data lives in AMReX managed memory so that the compute()
 * kernel works on both CPU and GPU builds.
 */

#ifndef EXTERNALGRID_INITIALDATA_HPP_
#define EXTERNALGRID_INITIALDATA_HPP_

#include "StateVariables.hpp"
#include <AMReX_Array4.H>
#include <AMReX_GpuMemory.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

class ExternalGridInitialData
{
  public:
    struct params_t
    {
        std::string gridinit_file;
        std::array<double, AMREX_SPACEDIM> grid_center{};
    };

    ExternalGridInitialData(const params_t &a_params, double a_dx)
        : m_dx(a_dx)
    {
        load_gridinit(a_params.gridinit_file);
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    compute(int i, int j, int k,
            amrex::Array4<amrex::Real> cell) const
    {
        const double px = (i + 0.5) * m_dx;
        const double py = (j + 0.5) * m_dx;
        const double pz = (k + 0.5) * m_dx;

        const double fi = (px - m_origin[0]) / m_file_dx[0] - 0.5;
        const double fj = (py - m_origin[1]) / m_file_dx[1] - 0.5;
        const double fk = (pz - m_origin[2]) / m_file_dx[2] - 0.5;

        const int i0 = amrex::max(0, amrex::min(int(fi), m_nx - 2));
        const int j0 = amrex::max(0, amrex::min(int(fj), m_ny - 2));
        const int k0 = amrex::max(0, amrex::min(int(fk), m_nz - 2));

        const double wx = amrex::max(0.0, amrex::min(1.0, fi - i0));
        const double wy = amrex::max(0.0, amrex::min(1.0, fj - j0));
        const double wz = amrex::max(0.0, amrex::min(1.0, fk - k0));

        const int nc = amrex::min(m_ncomp, static_cast<int>(NUM_VARS));
        for (int c = 0; c < nc; ++c)
        {
            const double v000 = m_ptr[flat(k0,     j0,     i0,     c)];
            const double v100 = m_ptr[flat(k0,     j0,     i0 + 1, c)];
            const double v010 = m_ptr[flat(k0,     j0 + 1, i0,     c)];
            const double v110 = m_ptr[flat(k0,     j0 + 1, i0 + 1, c)];
            const double v001 = m_ptr[flat(k0 + 1, j0,     i0,     c)];
            const double v101 = m_ptr[flat(k0 + 1, j0,     i0 + 1, c)];
            const double v011 = m_ptr[flat(k0 + 1, j0 + 1, i0,     c)];
            const double v111 = m_ptr[flat(k0 + 1, j0 + 1, i0 + 1, c)];

            const double val =
                v000 * (1 - wx) * (1 - wy) * (1 - wz) +
                v100 * wx       * (1 - wy) * (1 - wz) +
                v010 * (1 - wx) * wy       * (1 - wz) +
                v110 * wx       * wy       * (1 - wz) +
                v001 * (1 - wx) * (1 - wy) * wz +
                v101 * wx       * (1 - wy) * wz +
                v011 * (1 - wx) * wy       * wz +
                v111 * wx       * wy       * wz;

            cell(i, j, k, c) = static_cast<amrex::Real>(val);
        }

        cell(i, j, k, c_chi) =
            amrex::max(cell(i, j, k, c_chi), amrex::Real(1.0e-10));
        cell(i, j, k, c_lapse) =
            amrex::max(cell(i, j, k, c_lapse), amrex::Real(1.0e-10));
    }

  private:
    double m_dx{1.0};
    int m_nx{0}, m_ny{0}, m_nz{0}, m_ncomp{0};
    double m_file_dx[3]{1.0, 1.0, 1.0};
    double m_origin[3]{};

    // Managed memory: accessible from both host and device
    std::shared_ptr<double> m_managed;
    double *m_ptr{nullptr};

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE long
    flat(int kk, int jj, int ii, int c) const
    {
        return (static_cast<long>(kk) * m_ny + jj) * m_nx * m_ncomp +
               static_cast<long>(ii) * m_ncomp + c;
    }

    void load_gridinit(const std::string &path)
    {
        std::ifstream fin(path, std::ios::binary);
        if (!fin.is_open())
        {
            throw std::runtime_error(
                "ExternalGridInitialData: cannot open " + path);
        }

        std::string line;
        while (std::getline(fin, line))
        {
            if (line == "END_HEADER")
                break;

            std::istringstream iss(line);
            std::string key;
            iss >> key;
            if (key == "num_components")
                iss >> m_ncomp;
            else if (key == "nx_ny_nz")
                iss >> m_nx >> m_ny >> m_nz;
            else if (key == "dx")
            {
                iss >> m_file_dx[0];
                if (!(iss >> m_file_dx[1] >> m_file_dx[2]))
                {
                    // V1 format: single scalar dx
                    m_file_dx[1] = m_file_dx[2] = m_file_dx[0];
                }
            }
            else if (key == "origin")
                iss >> m_origin[0] >> m_origin[1] >> m_origin[2];
        }

        if (m_nx <= 0 || m_ny <= 0 || m_nz <= 0 || m_ncomp <= 0)
        {
            throw std::runtime_error(
                "ExternalGridInitialData: invalid header in " + path);
        }

        const long total =
            static_cast<long>(m_nz) * m_ny * m_nx * m_ncomp;

        // Allocate in managed memory for GPU accessibility
        double *raw = static_cast<double *>(
            amrex::The_Managed_Arena()->alloc(
                static_cast<std::size_t>(total) * sizeof(double)));
        m_managed = std::shared_ptr<double>(
            raw,
            [](double *p)
            { amrex::The_Managed_Arena()->free(p); });
        m_ptr = raw;

        fin.read(reinterpret_cast<char *>(m_ptr),
                 static_cast<std::streamsize>(total * sizeof(double)));

        if (!fin)
        {
            throw std::runtime_error(
                "ExternalGridInitialData: failed to read data from " +
                path);
        }

        amrex::Print() << "ExternalGridInitialData: loaded " << path
                        << " (" << m_nx << "x" << m_ny << "x" << m_nz
                        << ", " << m_ncomp << " components)\n";
    }
};

#endif /* EXTERNALGRID_INITIALDATA_HPP_ */
