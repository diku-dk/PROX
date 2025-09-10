#ifndef GEOMETRY_CLOSEST_POINT_ON_PLANE_H
#define GEOMETRY_CLOSEST_POINT_ON_PLANE_H

#include <types/geometry_plane.h>

#include <tiny_is_finite.h>
#include <tiny_is_number.h>

#include <cassert>

namespace geometry
{

template <typename T>
inline EigenVector3<T> closest_point_on_plane(const EigenVector3<T>& p,
                                              Plane<T> const& P)
{

    const EigenVector3<T>& n = P.n();
    T const & w = P.w();

    const EigenVector3<T> q = p - n * ( n).dot( (p - n*w) );

    assert( is_number(q(0)) || !"closest_point_on_line(): NaN encountered");
    assert( is_finite(q(0)) || !"closest_point_on_line(): Inf encountered");
    assert( is_number(q(1)) || !"closest_point_on_line(): NaN encountered");
    assert( is_finite(q(1)) || !"closest_point_on_line(): Inf encountered");
    assert( is_number(q(2)) || !"closest_point_on_line(): NaN encountered");
    assert( is_finite(q(2)) || !"closest_point_on_line(): Inf encountered");

    return q;
}

}// namespace geometry

// GEOMETRY_CLOSEST_POINT_ON_PLANE_H
#endif
