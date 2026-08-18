#ifndef AHFINDERSTATE_HPP_
#define AHFINDERSTATE_HPP_

#include <AMReX_IntegratorBase.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <concepts>
#include <cstddef>
#include <memory>
#include <vector>

// State of AHFinder (h and v for all particles) to be evolved by amrex
// timeintegrator
struct AHState
{
    std::vector<double> h;
    std::vector<double> v;

    AHState() = default;

    AHState(std::vector<double> a_h, std::vector<double> a_v)
        : h(std::move(a_h)), v(std::move(a_v))
    {
    }
};

namespace amrex
{

// Implementations of allocation, copying and SAXPY needed for time integration
template <class T>
    requires(std::same_as<T, AHState>)
struct IntegratorOps<T>
{
    static void CreateLike(amrex::Vector<std::unique_ptr<T>> &V, const T &Other,
                           bool /* Grow */ = false)
    {
        // Emplace a new, same-sized, zero-initialised state
        V.emplace_back(std::make_unique<T>());
        V.back()->h.assign(Other.h.size(), 0.0);
        V.back()->v.assign(Other.v.size(), 0.0);
    }

    static void Copy(T &Y, const T &Other)
    {
        // Copy the contents of Other into Y
        Y.h = Other.h;
        Y.v = Other.v;
    }

    static void Saxpy(T &Y, const amrex::Real a, const T &X)
    {
        // Calculate Y += a * X elementwise over both fields
        const std::size_t n = Y.h.size();
        for (std::size_t i = 0; i < n; ++i)
        {
            Y.h[i] += a * X.h[i];
            Y.v[i] += a * X.v[i];
        }
    }
};

} // namespace amrex

#endif /* AHFINDERSTATE_HPP_ */
