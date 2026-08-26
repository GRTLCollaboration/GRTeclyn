/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef INTERPOLATIONQUERYPARTICLE_HPP_
#define INTERPOLATIONQUERYPARTICLE_HPP_

// Other includes
#include "BCParity.hpp"
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
        int comp{};
        amrex::ParticleReal *out_data_ptr{};
        BCParity parity{
            BCParity::undefined}; // default parity is undefined, but should be
                                  // set for diagnostic variables
    };

    using comp_map_t = std::map<Derivative, std::vector<out_t>>;
    using iterator =
        typename std::map<Derivative, std::vector<out_t>>::const_iterator;

  private:
    template <int num_components> friend class ParticleInterpolator;

    size_t m_num_points;

    std::array<const amrex::ParticleReal *, AMREX_SPACEDIM> m_coords{};
    comp_map_t m_comps{};
    VariableType m_variable_type{}; // for a given InterpolationQueryParticle
                                    // the variable type must be the same!
    bool m_variable_type_set =
        false; // flag to check whether var type has been set

  public:
    InterpolationQueryParticle(int num_points) : m_num_points(num_points) {}

    // Returns the pointer that was passed to setCoords
    // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
    [[nodiscard]] const amrex::ParticleReal *coords(int dim) const

    {
        AMREX_ASSERT(dim >= 0 && dim < AMREX_SPACEDIM);
        return m_coords[dim];
    }

    InterpolationQueryParticle &setCoords(int dim,
                                          const amrex::ParticleReal *coords)
    {
        AMREX_ASSERT(dim < AMREX_SPACEDIM);
        this->m_coords[dim] = coords;
        return *this;
    }

    InterpolationQueryParticle &
    addComp(int comp, amrex::ParticleReal *out_ptr,
            VariableType variable_type = VariableType::state,
            BCParity parity            = BCParity::undefined,
            const Derivative &deriv    = Derivative::LOCAL)
    {
        AMREX_ASSERT(out_ptr != NULL || m_num_points == 0);

        if (!m_variable_type_set)
        {
            m_variable_type = variable_type; // first comp sets var type here
            m_variable_type_set = true;
        }
        else
        {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                m_variable_type == variable_type,
                "InterpolationQueryParticle::addComp(): VariableType mismatch "
                "across comps");
        }

        if (variable_type == VariableType::derived &&
            parity == BCParity::undefined)
        {
            amrex::Abort("InterpolationQueryParticle::addComp(): Parity has "
                         "not been defined! You should set it appropriately "
                         "for diagnostic variables!");
        }

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

        result->second.push_back(out_t{comp, out_ptr, parity});
        return *this;
    }

    InterpolationQueryParticle &clearComps()
    {
        m_comps.clear();
        m_variable_type_set = false; // reset the initialised var type flag
        return *this;
    }

    [[nodiscard]] VariableType getVariableType() const
    {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            m_variable_type_set, "InterpolationQueryParticle: VariableType "
                                 "accessed before any comps were added");
        return m_variable_type;
    }

    // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
    int numComps()
    {
        int accum = 0;

        for (auto &m_comp : m_comps)
        {
            accum += static_cast<int>(m_comp.second.size());
        }

        return accum;
    }

    [[nodiscard]] size_t numPoints() const { return m_num_points; }

    [[nodiscard]] iterator compsBegin() const { return m_comps.cbegin(); }

    [[nodiscard]] iterator compsEnd() const { return m_comps.cend(); }
};

#endif /* INTERPOLATIONQUERYPARTICLE_HPP_ */
