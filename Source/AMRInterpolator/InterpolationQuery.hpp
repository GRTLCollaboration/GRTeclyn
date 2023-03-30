/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INTERPOLATIONQUERY_HPP_
#define INTERPOLATIONQUERY_HPP_

// Other includes
#include "Derivative.hpp"
#include "VariableType.hpp"
#include <map>
#include <tuple>
#include <utility>
#include <vector>

class InterpolationQuery
{
  public:
    using out_t      = std::tuple<int, double *, VariableType>;
    using comp_map_t = std::map<Derivative, std::vector<out_t>>;
    using iterator =
        typename std::map<Derivative, std::vector<out_t>>::iterator;

  private:
    template <typename InterpAlgo> friend class AMRInterpolator;

    size_t m_num_points;
    std::vector<const double *> m_coords;
    comp_map_t m_comps;

  public:
    InterpolationQuery(int num_points)
        : m_num_points(num_points), m_coords(AMREX_SPACEDIM, nullptr)
    {
    }

    InterpolationQuery &setCoords(int dim, const double *coords)
    {
        AMREX_ASSERT(dim < AMREX_SPACEDIM);
        this->m_coords[dim] = coords;
        return *this;
    }

    InterpolationQuery &
    addComp(int comp, double *out_ptr,
            const Derivative &deriv    = Derivative::LOCAL,
            VariableType variable_type = VariableType::evolution)
    {
        AMREX_ASSERT(out_ptr != NULL || m_num_points == 0);

        auto result = m_comps.find(deriv);
        if (result == m_comps.end())
        {
            result = m_comps
                         .insert(std::pair<Derivative, std::vector<out_t>>(
                             deriv, std::vector<out_t>()))
                         .first;
        }

        result->second.emplace_back(comp, out_ptr, variable_type);
        return *this;
    }

    InterpolationQuery &clearComps()
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

#endif /* INTERPOLATIONQUERY_HPP_ */
