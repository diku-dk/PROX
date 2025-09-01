#ifndef GEOMETRY_BARYCENTRIC_INSIDE_H
#define GEOMETRY_BARYCENTRIC_INSIDE_H

#include <tiny_value_traits.h>
#include <tiny_is_number.h>
#include <tiny_is_finite.h>

#include <cmath>
#include <cassert>

namespace geometry
{

  /**
   * Test of barycentric coordinates of an edge corresponds to a point inside the edge.
   */
  template<typename T>
  inline bool barycentric_inside(T const & w1, T const & w2)
  {
    assert( is_number(w1) || !"barycentric_inside(): NaN encountered");
    assert( is_number(w2) || !"barycentric_inside(): NaN encountered");
    assert( is_finite(w1) || !"barycentric_inside(): Inf encountered");
    assert( is_finite(w2) || !"barycentric_inside(): Inf encountered");

    if(w1<0)
      return false;
    if(w2<0)
      return false;

    if(w1>1)
      return false;
    if(w2>1)
      return false;

    if((w1+w2)>1)
      return false;

    if((w1+w2)<0)
      return false;

    return true;
  }

  /**
   * Test of barycentric coordinates of a triangle corresponds to a point inside the triangle.
   */
  template<typename T>
  inline bool barycentric_inside(T const & w1, T const & w2, T const & w3)
  {
    assert( is_number(w1) || !"barycentric_inside(): NaN encountered");
    assert( is_number(w2) || !"barycentric_inside(): NaN encountered");
    assert( is_number(w3) || !"barycentric_inside(): NaN encountered");
    assert( is_finite(w1) || !"barycentric_inside(): Inf encountered");
    assert( is_finite(w2) || !"barycentric_inside(): Inf encountered");
    assert( is_finite(w3) || !"barycentric_inside(): Inf encountered");


    if(w1<0)
      return false;
    if(w2<0)
      return false;
    if(w3<0)
      return false;

    if(w1>1)
      return false;
    if(w2>1)
      return false;
    if(w3>1)
      return false;

    if((w1+w2+w3)>1)
      return false;

    if((w1+w2+w3)<0)
      return false;


    return true;
  }

  /**
   * Test of barycentric coordinates of a tetrahedron corresponds to a point inside the tetrahedron.
   */
  template<typename T>
  inline bool barycentric_inside(T const & w1, T const & w2, T const & w3, T const & w4)
  {
    typedef typename tiny::ValueTraits<T> VT;

    assert( is_number(w1) || !"barycentric_inside(): NaN encountered");
    assert( is_number(w2) || !"barycentric_inside(): NaN encountered");
    assert( is_number(w3) || !"barycentric_inside(): NaN encountered");
    assert( is_number(w4) || !"barycentric_inside(): NaN encountered");
    assert( is_finite(w1) || !"barycentric_inside(): Inf encountered");
    assert( is_finite(w2) || !"barycentric_inside(): Inf encountered");
    assert( is_finite(w3) || !"barycentric_inside(): Inf encountered");
    assert( is_finite(w4) || !"barycentric_inside(): Inf encountered");

    if(w1<0)
      return false;
    if(w2<0)
      return false;
    if(w3<0)
      return false;
    if(w4<0)
      return false;

    if(w1>1)
      return false;
    if(w2>1)
      return false;
    if(w3>1)
      return false;
    if(w4>1)
      return false;

    if((w1+w2+w3+w4)>1)
      return false;

    if((w1+w2+w3+w4)<0)
      return false;

    return true;
  }

}//namespace geometry

// GEOMETRY_BARYCENTRIC_INSIDE_H
#endif
