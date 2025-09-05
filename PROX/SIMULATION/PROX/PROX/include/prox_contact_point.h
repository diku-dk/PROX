#ifndef PROX_CONTACT_POINT_H
#define PROX_CONTACT_POINT_H

#include "contacts/contact_point.h"
#include <prox_rigid_body.h>

#include <tiny_precision.h>   // needed for working precision
#include <tiny_math_types.h>

#include <cmath> // needed for std::fabs
#include <cassert> // needed for assert

namespace prox
{

template <typename M> using ContactPoint = geometry::ContactPoint<RigidBody<M>, typename RigidBody<M>::T>;

}// namespace prox

// PROX_CONTACT_POINT_H
#endif
