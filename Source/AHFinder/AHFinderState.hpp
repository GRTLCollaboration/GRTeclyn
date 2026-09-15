#ifndef AHFINDERSTATE_HPP_
#define AHFINDERSTATE_HPP_

#include <utility>
#include <vector>

// State of AHFinder: the surface radius h for all grid points. The implicit
// pseudo-transient continuation (PTC) scheme integrates dh/dt = -Theta(h)
// directly, so the state is just h.
struct AHState
{
    std::vector<double> h;

    AHState() = default;

    explicit AHState(std::vector<double> a_h) : h(std::move(a_h)) {}
};

#endif /* AHFINDERSTATE_HPP_ */
