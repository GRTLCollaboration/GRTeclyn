/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef INTERPOLATIONQUERYPARTICLE_HPP_
#define INTERPOLATIONQUERYPARTICLE_HPP_

// Other includes
#include "Derivative.hpp"
#include "VariableType.hpp"
#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

class InterpolationQueryParticle
{
  public:
    struct out_t
    {
        int comp;
        amrex::ParticleReal *out_data_ptr;
        VariableType variable_type; // state or derived?
    };

    using comp_map_t = std::map<Derivative, std::vector<out_t>>;
    using iterator =
        typename std::map<Derivative, std::vector<out_t>>::iterator;

  private:
    template <int num_components> friend class ParticleInterpolator;

    size_t m_num_points;
    std::array<const double *, AMREX_SPACEDIM> m_coords;
    comp_map_t m_comps;

  public:
    InterpolationQueryParticle(int num_points)
        : m_num_points(num_points), m_coords{}
    {
    }

    // Returns the pointer that was passed to setCoords
    inline const double *coords(int dim) const
    {
        AMREX_ASSERT(dim >= 0 && dim < AMREX_SPACEDIM);
        return m_coords[dim];
    }

    InterpolationQueryParticle &setCoords(int dim, const double *coords)
    {
        AMREX_ASSERT(dim < AMREX_SPACEDIM);
        this->m_coords[dim] = coords;
        return *this;
    }

    InterpolationQueryParticle &
    addComp(int comp, double *out_ptr,
            const Derivative &deriv    = Derivative::LOCAL,
            VariableType variable_type = VariableType::state)
    {
        AMREX_ASSERT(out_ptr != NULL || m_num_points == 0);

        // for now we do not allow derivatives
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
        {
            if (deriv[dir] != 0)
            {
                amrex::Abort(
                    "InterpolationQueryParticle::addComp(): "
                    "Derivative interpolation is not yet implemented :/ !");
            }
        }

        auto result = m_comps.find(deriv);
        if (result == m_comps.end())
        {
            result = m_comps
                         .insert(std::pair<Derivative, std::vector<out_t>>(
                             deriv, std::vector<out_t>()))
                         .first;
        }

        result->second.push_back(out_t{comp, out_ptr, variable_type});
        return *this;
    }

    InterpolationQueryParticle &clearComps()
    {
        m_comps.clear();
        return *this;
    }

    inline int numComps()
    {
        int accum = 0;

        for (auto &m_comp : m_comps)
        {
            accum += static_cast<int>(m_comp.second.size());
        }

        return accum;
    }

    [[nodiscard]] inline size_t numPoints() const { return m_num_points; }

    inline iterator compsBegin() { return m_comps.begin(); }

    inline iterator compsEnd() { return m_comps.end(); }
};

#endif /* INTERPOLATIONQUERYPARTICLE_HPP_ */
