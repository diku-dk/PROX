#ifndef PROX_GET_PRE_STABILIZATION_VECTOR_H
#define PROX_GET_PRE_STABILIZATION_VECTOR_H

#include "eigenhelperfunctions.h"
#include <steppers/prox_stepper_params.h>

#include <tiny_is_number.h>

#include <util_log.h>

#include <cassert>
#include <cmath>

/**
   *
   * @param drift_reduction     The procentage of numerical drift reduction that should
   *                            be done. This parameter must be in the range [0..1].
   * @param time_step           The time step-length. This paramter is needed to convert
   *                            drift-correction from displacement level to velocity level.
   *
   * @param w                  Current relative contact point velocites, w = J*u
   *
   */
namespace prox
{
template <typename T, typename Iterator>
inline void get_pre_stabilization_vector_eigen(Iterator begin, Iterator end,
                                               const StepperParams<T>& params,
                                               T timestep,
                                               const Eigen::VectorX<T>& w,
                                               Eigen::VectorX<T>& g, size_t K)
{
    assert(timestep > 0
           && "get_pre_stabilization_vector(): time_step should be positive");

    util::Log logging;

    g.resize(K * 4);  // Resize to (K * 4) x 1
    g.setZero();      // Set all elements to 0

    logging << "get_pre_stabilization_vector(): pre stabilizaiton = "
            << params.pre_stabilization() << util::Log::newline();

    if (!params.pre_stabilization()) return;

    T const& reduction = params.gap_reduction();
    T const& min_gap = params.min_gap();
    T const& max_gap = params.max_gap();

    assert(reduction >= 0
           && "get_pre_stabilization_vector(): gap reduction parameter should "
              "be positive");
    assert(reduction <= 1
           && "get_pre_stabilization_vector(): gap reduction parameter should "
              "be less than or "
              "equal to one");
    assert(min_gap >= 0
           && "get_pre_stabilization_vector(): min gap correction should be "
              "non negative");
    assert(max_gap > 0
           && "get_pre_stabilization_vector(): max gap correction should be "
              "positive");

    T const k = reduction / timestep;
    T const limit = -max_gap / timestep;
    T const yield = -min_gap;

    size_t index = 0;
    for (auto contact = begin; contact != end; ++contact, ++index)
    {
        // Calculate indices for this contact
        size_t w_idx = index * 4;
        size_t g_idx = index * 4;

        // Get the normal velocity component
        T const& v_n = w(w_idx);
        bool const add_correction = contact->depth <= yield && v_n <= 0;

        // Set the first component of this contact's data
        g(g_idx) = add_correction ? std::max(limit, k * contact->depth) : 0;

        // The other three components remain 0 (already set by setZero())

        assert(is_number(g(g_idx))
               && "get_pre_stabilization_vector(): g(g_idx) is not a number");
    }
}
} //namespace prox

// PROX_GET_PRE_STABILIZATION_VECTOR_H
#endif
