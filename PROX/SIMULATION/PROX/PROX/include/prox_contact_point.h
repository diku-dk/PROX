#ifndef PROX_CONTACT_POINT_H
#define PROX_CONTACT_POINT_H

#include "contacts/contact_point.h"
#include <prox_rigid_body.h>

#include <eigenhelperall.h>

#include <cmath> // needed for std::fabs
#include <cassert> // needed for assert

namespace prox
{

template <typename T>
requires(std::is_floating_point_v<T>)
using ContactPoint = geometry::ContactPoint<RigidBody<T>, T>;

}// namespace prox

// PROX_CONTACT_POINT_H
#endif
