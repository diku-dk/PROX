#ifndef HYPER_CONTACT_POINT_H
#define HYPER_CONTACT_POINT_H

#include "contacts/contact_point.h"
#include <hyper_body.h>

namespace hyper
{

template <typename MT> using ContactPoint = geometry::ContactPoint<Body<MT>, typename Body<MT>::T>;

}// namespace hyper

// HYPER_CONTACT_POINT_H
#endif
