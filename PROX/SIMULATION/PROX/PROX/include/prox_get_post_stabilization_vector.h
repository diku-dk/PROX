#ifndef PROX_GET_POST_STABILIZATION_VECTOR_H
#define PROX_GET_POST_STABILIZATION_VECTOR_H

#include "eigenhelperfunctions.h"
#include <steppers/prox_stepper_params.h>

#include <tiny_is_number.h>

#include <util_log.h>

#include <cassert>
#include <cmath>

namespace prox
{

template <typename T, typename Iterator>
inline void get_post_stabilization_vector_eigen(Iterator begin, Iterator end,
                                                const StepperParams<T>& params,
                                                Eigen::VectorX<T>& g, size_t K)
{
    util::Log logging;

    g.resize(K * 4);  // Resize to (K*4)x1 vector

    logging << "get_post_stabilization_vector(): post stabilization = "
            << params.post_stabilization() << util::Log::newline();

    if (!params.post_stabilization())
    {
        g.setZero();  // Set all elements to zero if no post stabilization
        return;
    }

    T const& reduction = params.gap_reduction();
    T const& max_gap = params.max_gap();

    assert(
        reduction >= 0
        || !"get_post_stabilization_vector(): reduction must be non-negative");
    assert(reduction <= 1
           || !"get_post_stabilization_vector(): reduction must less than or "
               "equal one");
    assert(
        max_gap > 0
        || !"get_post_stabilization_vector(): max reduction must be positive");

    // Initialize all elements to zero
    g.setZero();

    size_t index = 0u;

    for (auto contact = begin; contact != end; ++contact, ++index)
    {
        // Calculate the position in the flat vector
        size_t pos = index * 4;

        g(pos) = std::max(-max_gap, std::min<T>(reduction * contact->depth, 0));

        assert(is_number(g(pos))
               || !"get_post_stabilization_vector(): g(pos) is not a number");
    }
}
} //namespace prox

// PROX_GET_POST_STABILIZATION_VECTOR_H
#endif
