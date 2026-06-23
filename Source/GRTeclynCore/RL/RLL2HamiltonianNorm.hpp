#ifndef RL_L2_HAMILTONIAN_NORM_HPP_
#define RL_L2_HAMILTONIAN_NORM_HPP_

#include "RLMatterPumpParams.hpp"

//! Shared L2 Hamiltonian cache publisher (reducer lives in example code).
inline void publish_cached_L2_Ham(double l2_ham)
{
    RLRuntime::publish_cached_L2_Ham(l2_ham);
}

#endif /* RL_L2_HAMILTONIAN_NORM_HPP_ */
