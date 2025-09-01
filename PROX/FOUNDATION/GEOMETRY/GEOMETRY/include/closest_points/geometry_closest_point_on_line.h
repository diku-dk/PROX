#ifndef GEOMETRY_CLOSEST_POINT_ON_LINE_H
#define GEOMETRY_CLOSEST_POINT_ON_LINE_H

#include <types/geometry_line.h>

#include <tiny_is_finite.h>
#include <tiny_is_number.h>



#include <cassert>


namespace geometry
{

  template<typename T>
  inline EigenVector3<T> closest_point_on_line(const EigenVector3<T>& p, Line<T> const & L)
  {
    const EigenVector3<T>& o = L.point();
    const EigenVector3<T>& d = L.direction();

    const  EigenVector3<T> q = o + d * ( d).dot( (p-o) );

    assert( is_number(q(0)) || !"closest_point_on_line(): NaN encountered");
    assert( is_finite(q(0)) || !"closest_point_on_line(): Inf encountered");
    assert( is_number(q(1)) || !"closest_point_on_line(): NaN encountered");
    assert( is_finite(q(1)) || !"closest_point_on_line(): Inf encountered");
    assert( is_number(q(2)) || !"closest_point_on_line(): NaN encountered");
    assert( is_finite(q(2)) || !"closest_point_on_line(): Inf encountered");

    return q;
  }

}// namespace geometry

// GEOMETRY_CLOSEST_POINT_ON_LINE_H
#endif
