#ifndef PROX_UPDATE_INERTIA_TENSOR_H
#define PROX_UPDATE_INERTIA_TENSOR_H

#include "tiny_math_types.h"

#include <type_traits>

namespace prox
{
namespace detail
{

/**
     * Inertia Update Method.
     *
     * Computes, I_world = R I_body R^T. This method have been optimized to
     * avoid having to transpose the orientation matrix and exploit
     * symmetry of the inertia tensor.
     *
     * @note            Further optimization may be possible by exploiting common sub-terms.
     *
     * @param R         The orientation of the body frame with repect the world frame as a rotation matrix.
     * @param I_body    An inertia tensor with respect to the body frame
     * @param I_world   Upon return this argument holds the inertia tensor with respect to the orientation of the world frame.
     *
     * @tparam M        The 3-by-3 matrix type to use for this update.
     */
template <typename T>
requires std::is_floating_point_v<T>
inline void update_inertia_tensor(EigenMatrix3<T>& R, const EigenMatrix3<T>& I_body,
                                  EigenMatrix3<T>& I_world, const bool = false)
{

    //--- 2004-08-30 Kenny: The formulas below was generated using Matlab
    //---
    //--- syms R00 R01 R02  R10 R11 R12 R20 R21 R22 real;
    //--- syms I00 I01 I02  I10 I11 I12 I20 I21 I22 real;
    //--- I = [ I00, I01 , I02; I01, I11, I12; I02, I12, I22]
    //--- R = [ R00, R01 , R02; R10, R11, R12; R20, R21, R22]
    //--- W = R*I*R'
    //--- simplify(W)
    //---

    I_world(0, 0) = R(0, 0) * R(0, 0) * I_body(0, 0) + 2 * R(0, 0) * R(0, 1) * I_body(0, 1)
                  + 2 * R(0, 0) * R(0, 2) * I_body(0, 2) + R(0, 1) * R(0, 1) * I_body(1, 1)
                  + 2 * R(0, 1) * R(0, 2) * I_body(1, 2) + R(0, 2) * R(0, 2) * I_body(2, 2);
    I_world(1, 1) = R(1, 0) * R(1, 0) * I_body(0, 0) + 2 * R(1, 0) * R(1, 1) * I_body(0, 1)
                  + 2 * R(1, 0) * R(1, 2) * I_body(0, 2) + R(1, 1) * R(1, 1) * I_body(1, 1)
                  + 2 * R(1, 1) * R(1, 2) * I_body(1, 2) + R(1, 2) * R(1, 2) * I_body(2, 2);
    I_world(2, 2) = R(2, 0) * R(2, 0) * I_body(0, 0) + 2 * R(2, 0) * R(2, 1) * I_body(0, 1)
                  + 2 * R(2, 0) * R(2, 2) * I_body(0, 2) + R(2, 1) * R(2, 1) * I_body(1, 1)
                  + 2 * R(2, 1) * R(2, 2) * I_body(1, 2) + R(2, 2) * R(2, 2) * I_body(2, 2);
    I_world(1, 0) = I_world(0, 1)
        = R(1, 0) * R(0, 0) * I_body(0, 0) + R(1, 0) * R(0, 1) * I_body(0, 1)
        + R(1, 0) * R(0, 2) * I_body(0, 2) + R(1, 1) * R(0, 0) * I_body(0, 1)
        + R(1, 1) * R(0, 1) * I_body(1, 1) + R(1, 1) * R(0, 2) * I_body(1, 2)
        + R(1, 2) * R(0, 0) * I_body(0, 2) + R(1, 2) * R(0, 1) * I_body(1, 2)
        + R(1, 2) * R(0, 2) * I_body(2, 2);
    I_world(2, 0) = I_world(0, 2)
        = R(2, 0) * R(0, 0) * I_body(0, 0) + R(2, 0) * R(0, 1) * I_body(0, 1)
        + R(2, 0) * R(0, 2) * I_body(0, 2) + R(2, 1) * R(0, 0) * I_body(0, 1)
        + R(2, 1) * R(0, 1) * I_body(1, 1) + R(2, 1) * R(0, 2) * I_body(1, 2)
        + R(2, 2) * R(0, 0) * I_body(0, 2) + R(2, 2) * R(0, 1) * I_body(1, 2)
        + R(2, 2) * R(0, 2) * I_body(2, 2);
    I_world(2, 1) = I_world(1, 2)
        = R(1, 0) * R(2, 0) * I_body(0, 0) + R(1, 0) * R(2, 1) * I_body(0, 1)
        + R(1, 0) * R(2, 2) * I_body(0, 2) + R(1, 1) * R(2, 0) * I_body(0, 1)
        + R(1, 1) * R(2, 1) * I_body(1, 1) + R(1, 1) * R(2, 2) * I_body(1, 2)
        + R(1, 2) * R(2, 0) * I_body(0, 2) + R(1, 2) * R(2, 1) * I_body(1, 2)
        + R(1, 2) * R(2, 2) * I_body(2, 2);
}

} // namespace detail
} // namespace prox

// PROX_UPDATE_INERTIA_TENSOR_H
#endif
