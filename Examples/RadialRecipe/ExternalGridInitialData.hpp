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
 *
 * COMPONENT MATCHING IS BY NAME.
 * This class used to copy file component c into state variable c and never read
 * the header's `component_names` at all, truncating at min(ncomp, NUM_VARS)
 * without a word.  The Python side of the same handoff indexes the file *by
 * name* (`comp_names.index("lapse")`).  Two codes, two different notions of what
 * identifies a component, and nothing comparing them: reorder or insert one
 * state variable and every field lands one slot off, silently, producing a
 * plausible-looking spacetime that solves nothing.  See DebugPreGPU.md.
 *
 * We now resolve every component through StateVariables::names and refuse to
 * load a file we cannot map exactly.  Files predating `component_names` still
 * load positionally, but say so loudly.
 */

#ifndef EXTERNALGRID_INITIALDATA_HPP_
#define EXTERNALGRID_INITIALDATA_HPP_

#include "StateVariables.hpp"
#include <AMReX_Array4.H>
#include <AMReX_GpuMemory.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

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

        // Iterate over DESTINATION slots: m_comp_of_var[v] is the file column
        // that feeds state variable v, resolved by name at load time, or -1 if
        // the file does not carry that variable (left at its pre-set value).
        for (int v = 0; v < static_cast<int>(NUM_VARS); ++v)
        {
            const int c = m_comp_of_var[v];
            if (c < 0)
                continue;

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

            cell(i, j, k, v) = static_cast<amrex::Real>(val);
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

    //! For each GRTeclyn state variable, the file column that feeds it, or -1
    //! when the file does not carry it.  Built by name in load_gridinit(); a
    //! trivially-copyable plain array so it travels into device lambdas.
    int m_comp_of_var[NUM_VARS];

    // Managed memory: accessible from both host and device
    std::shared_ptr<double> m_managed;
    double *m_ptr{nullptr};

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE long
    flat(int kk, int jj, int ii, int c) const
    {
        return (static_cast<long>(kk) * m_ny + jj) * m_nx * m_ncomp +
               static_cast<long>(ii) * m_ncomp + c;
    }

    //! Names the GRTresna bridge writes that differ from StateVariables::names
    //! for the same slot.  StateVariables.hpp:12-13: the single-complex-scalar
    //! model reuses the first lump slots, and the writer emits the enum
    //! spelling (c_phi2 / c_Pi2) while names[] carries the lump spelling.  Both
    //! are correct for the same component; neither is a mismatch.
    static std::string canonical_component_name(const std::string &name)
    {
        if (name == "phi2")
            return "phi_lump0";
        if (name == "Pi2")
            return "Pi_lump0";
        return name;
    }

    //! Resolve every file column to a state-variable slot BY NAME.
    void build_component_map(const std::vector<std::string> &comp_names,
                             const std::string &path)
    {
        for (int v = 0; v < static_cast<int>(NUM_VARS); ++v)
            m_comp_of_var[v] = -1;

        if (comp_names.empty())
        {
            // Pre-V2 file: no names to match on, so fall back to the positional
            // contract this class used to apply unconditionally.  Say so --
            // nothing downstream can tell a correct load from a shifted one.
            const int nc = std::min(m_ncomp, static_cast<int>(NUM_VARS));
            for (int c = 0; c < nc; ++c)
                m_comp_of_var[c] = c;
            amrex::Print()
                << "WARNING: ExternalGridInitialData: " << path
                << " carries no component_names; matching by POSITION. A "
                   "reordered or inserted state variable would be loaded into "
                   "the wrong slot with no error.\n";
            return;
        }

        if (static_cast<int>(comp_names.size()) != m_ncomp)
        {
            throw std::runtime_error(
                "ExternalGridInitialData: " + path + " declares num_components " +
                std::to_string(m_ncomp) + " but lists " +
                std::to_string(comp_names.size()) +
                " component_names. The payload stride is whichever one is "
                "wrong, so every row after the first would be misaligned.");
        }

        int n_mapped   = 0;
        int n_reordered = 0;
        for (int c = 0; c < m_ncomp; ++c)
        {
            const std::string name = canonical_component_name(comp_names[c]);

            int slot = -1;
            for (int v = 0; v < static_cast<int>(NUM_VARS); ++v)
            {
                if (StateVariables::names[v] == name)
                {
                    slot = v;
                    break;
                }
            }
            if (slot < 0)
            {
                throw std::runtime_error(
                    "ExternalGridInitialData: " + path + " component " +
                    std::to_string(c) + " is named '" + comp_names[c] +
                    "', which is not a GRTeclyn state variable. Refusing to "
                    "guess which slot it belongs in.");
            }
            if (m_comp_of_var[slot] >= 0)
            {
                throw std::runtime_error(
                    "ExternalGridInitialData: " + path + " lists '" +
                    comp_names[c] + "' twice (components " +
                    std::to_string(m_comp_of_var[slot]) + " and " +
                    std::to_string(c) + ").");
            }
            m_comp_of_var[slot] = c;
            ++n_mapped;
            if (slot != c)
                ++n_reordered;
        }

        amrex::Print() << "ExternalGridInitialData: mapped " << n_mapped
                       << " of " << static_cast<int>(NUM_VARS)
                       << " state variables by name";
        if (n_reordered > 0)
        {
            // Worth shouting about: under the old positional loader every one
            // of these landed in the wrong variable without complaint.
            amrex::Print() << " (" << n_reordered
                           << " would have been MISPLACED by position)";
        }
        amrex::Print() << "\n";
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
        std::vector<std::string> comp_names;
        while (std::getline(fin, line))
        {
            if (line == "END_HEADER")
                break;

            std::istringstream iss(line);
            std::string key;
            iss >> key;
            if (key == "num_components")
                iss >> m_ncomp;
            else if (key == "component_names")
            {
                std::string name;
                while (iss >> name)
                    comp_names.push_back(name);
            }
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

        build_component_map(comp_names, path);

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
