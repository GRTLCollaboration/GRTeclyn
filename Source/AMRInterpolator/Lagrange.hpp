/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef LAGRANGE_HPP_
#define LAGRANGE_HPP_

#include "InterpSource.hpp"
#include <utility>

template <int Order> class Lagrange
{
    const InterpSource &m_source;
    bool m_verbosity;

    struct Stencil;

    struct Stencil
    {
        int m_width;
        int m_deriv;
        double m_dx;
        double m_point_offset;

        std::vector<double> m_weights;

        Stencil(int width, int deriv, double dx, double point_offset);
        inline bool operator==(const Stencil &rhs) const;
        [[nodiscard]] inline bool isSameAs(int width, int deriv, double dx,
                                           double point_offset) const;

        inline const double &operator[](unsigned int i) const;
    };

    using stencil_collection_t = std::vector<Stencil>;
    stencil_collection_t m_memoized_stencils;

    Stencil getStencil(int width, int deriv, double dx, double point_offset);

    // Helper function to generate tensor product weights
    // Argument 'dim' is used for recursion over dimensions.
    std::pair<std::vector<amrex::IntVect>, std::vector<double>>
    generateStencil(const std::array<int, AMREX_SPACEDIM> &deriv,
                    const std::array<double, AMREX_SPACEDIM> &dx,
                    const std::array<double, AMREX_SPACEDIM> &evalCoord,
                    const amrex::IntVect &nearest,
                    int dim = AMREX_SPACEDIM - 1);

    std::vector<amrex::IntVect> m_interp_points;
    std::vector<double> m_interp_weights;

    // We are adding 216+ numbers at roughly the same magnitudes but alternating
    // signs. Let's keep track of positive and negative terms separately to make
    // sure we don't run into trouble.
    std::multiset<double> m_interp_neg;
    std::multiset<double> m_interp_pos;

  public:
    Lagrange(const InterpSource &source, bool verbosity = false);

    void setup(const std::array<int, AMREX_SPACEDIM> &deriv,
               const std::array<double, AMREX_SPACEDIM> &dx,
               const std::array<double, AMREX_SPACEDIM> &evalCoord,
               const amrex::IntVect &nearest);
    double interpData(const amrex::FArrayBox &fab, int comp);

    const static std::string TAG;
};

#include "Lagrange.impl.hpp"

#endif /* LAGRANGE_HPP_ */
