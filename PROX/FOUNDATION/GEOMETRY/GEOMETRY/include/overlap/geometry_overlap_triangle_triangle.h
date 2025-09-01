#ifndef GEOMETRY_OVERLAP_TRIANGLE_TRIANGLE_H
#define GEOMETRY_OVERLAP_TRIANGLE_TRIANGLE_H

#include <types/geometry_triangle.h>
#include <types/geometry_plane.h>

#include <tiny_vector_functions.h>

namespace geometry
{
  namespace detail
  {

    template< typename T>
    inline bool SAT_line( const  EigenVector3<T>& a0
                         , const EigenVector3<T>& a1
                         , const EigenVector3<T>& b0
                         , const EigenVector3<T>& b1
                         , const EigenVector3<T>& n
                         )
    {
      using std::max;
      using std::min;

      T const da0 = dot(n, a0);
      T const da1 = dot(n, a1);
      T const db0 = dot(n, b0);
      T const db1 = dot(n, b1);

      T const a_min = min(da0, da1);
      T const a_max = max(da0, da1);
      T const b_min = min(db0, db1);
      T const b_max = max(db0, db1);

      if( (a_max < b_min) || (a_min > b_max) )
        return true;

      return false;
    }

  } // end of namespace detail

  template< typename T>
  inline bool overlap_triangle_triangle(
                                        Triangle<T> const & A
                                        , Triangle<T> const & B
                                        )
  {

    const EigenVector3<T>& A0 = A.point(0);
    const EigenVector3<T>& A1 = A.point(1);
    const EigenVector3<T>& A2 = A.point(2);

    const EigenVector3<T>& B0 = B.point(0);
    const EigenVector3<T>& B1 = B.point(1);
    const EigenVector3<T>& B2 = B.point(2);

    // All vertices of B is on the same side of A
    {
      Plane<T>    const plane    = make_plane(A);

      T const b0 = get_signed_distance(B0, plane);
      T const b1 = get_signed_distance(B1, plane);
      T const b2 = get_signed_distance(B2, plane);

      if ( b0 > 0 && b1 > 0 && b2 > 0  )
        return false;

      if ( b0 < 0 && b1 < 0 && b2 < 0  )
        return false;
    }

    // All vertices of A is on the same side of B
    {
      Plane<T>    const plane    = make_plane(B);

      T const a0 = get_signed_distance(A0, plane);
      T const a1 = get_signed_distance(A1, plane);
      T const a2 = get_signed_distance(A2, plane);

      if ( a0 > 0 && a1 > 0 && a2 > 0  )
        return false;

      if ( a0 < 0 && a1 < 0 && a2 < 0  )
        return false;
    }

    // None of the face normals was a separating axis, now we test the edges

    // let edge0 be vertices 0 and 1
    // let edge1 be vertices 0 and 2
    // let edge2 be vertices 1 and 2

    // Ae0 v Be0
    const EigenVector3<T> ae0be0 = cross(A1-A0, B1-B0);
    if (detail::SAT_line(A1, A2, B1, B2, ae0be0))
      return false;

    // Ae0 v Be1
    const EigenVector3<T> ae0be1 = cross(A1-A0, B2-B0);
    if (detail::SAT_line(A1, A2, B0, B1, ae0be1))
      return false;

    // Ae0 v Be2
    const EigenVector3<T> ae0be2 = cross(A1-A0, B2-B1);
    if (detail::SAT_line(A1, A2, B0, B1, ae0be2))
      return false;

    // Ae1 v Be0
    const EigenVector3<T> ae1be0 = cross(A2-A0, B1-B0);
    if (detail::SAT_line(A1, A2, B1, B2, ae1be0))
      return false;

    // Ae1 v Be1
    const EigenVector3<T> ae1be1 = cross(A2-A0, B2-B0);
    if (detail::SAT_line(A1, A2, B1, B2, ae1be1))
      return false;

    // Ae1 v Be2
    const EigenVector3<T> ae1be2 = cross(A2-A0, B2-B1);
    if (detail::SAT_line(A1, A2, B0, B1, ae1be2))
      return false;

    // Ae2 v Be0
    const EigenVector3<T> ae2be0 = cross(A2-A1, B1-B0);
    if (detail::SAT_line(A0, A2, B1, B2, ae2be0))
      return false;

    // Ae2 v Be1
    const EigenVector3<T> ae2be1 = cross(A2-A1, B2-B0);
    if (detail::SAT_line(A0, A2, B1, B2, ae2be1))
      return false;

    // Ae2 v Be2
    const EigenVector3<T> ae2be2 = cross(A2-A1, B2-B1);
    if (detail::SAT_line(A0, A2, B0, B2, ae2be2))
      return false;

    return true;
  }

}// namespace geometry

// GEOMETRY_OVERLAP_TRIANGLE_TRIANGLE_H
#endif
