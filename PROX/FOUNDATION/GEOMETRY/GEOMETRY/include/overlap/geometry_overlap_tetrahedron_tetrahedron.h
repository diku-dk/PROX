#ifndef GEOMETRY_OVERLAP_TETRAHEDRON_TETRAHEDRON_H
#define GEOMETRY_OVERLAP_TETRAHEDRON_TETRAHEDRON_H

#include <types/geometry_tetrahedron.h>
#include <types/geometry_plane.h>

#include <tiny_vector_functions.h>

namespace geometry
{

    namespace detail
    {

    template< typename T>
    inline bool SAT_line( const  EigenVector3<T>& t0_p0
                         , const EigenVector3<T>& t0_p1
                         , const EigenVector3<T>& t0_p2
                         , const EigenVector3<T>& t1_p0
                         , const EigenVector3<T>& t1_p1
                         , const EigenVector3<T>& t1_p2
                         , const EigenVector3<T>& n
                         )
    {
        using std::max;
        using std::min;

        T const t0_d0 = dot(n, t0_p0);
        T const t0_d1 = dot(n, t0_p1);
        T const t0_d2 = dot(n, t0_p2);

        T const t1_d0 = dot(n, t1_p0);
        T const t1_d1 = dot(n, t1_p1);
        T const t1_d2 = dot(n, t1_p2);

        T const t0_min = min( min(t0_d0, t0_d1), t0_d2);
        T const t1_min = min( min(t1_d0, t1_d1), t1_d2);
        T const t0_max = max( max(t0_d0, t0_d1), t0_d2);
        T const t1_max = max( max(t1_d0, t1_d1), t1_d2);

        if( (t0_max < t1_min) || (t0_min > t1_max) )
            return true;

        return false;
    }

    } // end of namespace detail

  template< typename T>
  inline bool overlap_tetrahedron_tetrahedron(
                                              TetrahedronEigen<T> const & A
                                              , TetrahedronEigen<T> const & B
                                              )
  {

    assert( get_signed_volume(A) > 0 || !"overlap_tetrahedron_tetrahedron(): A was left oriented ");
    assert( get_signed_volume(B) > 0 || !"overlap_tetrahedron_tetrahedron(): B was left oriented ");

    const EigenVector3<T>& A0 = (A.point(0));
    const EigenVector3<T>& A1 = (A.point(1));
    const EigenVector3<T>& A2 = (A.point(2));
    const EigenVector3<T>& A3 = (A.point(3));

    const EigenVector3<T>& B0 = (B.point(0));
    const EigenVector3<T>& B1 = (B.point(1));
    const EigenVector3<T>& B2 = (B.point(2));
    const EigenVector3<T>& B3 = (B.point(3));

    for (unsigned int i = 0u; i < 4u; ++i)
    {
      Triangle<T> const triangle = get_opposite_face( i, A );
      Plane<T>    const plane    = make_plane(triangle);

      T const b0 = get_signed_distance(B0, plane);
      T const b1 = get_signed_distance(B1, plane);
      T const b2 = get_signed_distance(B2, plane);
      T const b3 = get_signed_distance(B3, plane);

      if ( b0 > 0 && b1 > 0 && b2 > 0 && b3 > 0 )
        return false;
    }

    for (unsigned int i = 0u; i < 4u; ++i)
    {
      Triangle<T> const triangle = get_opposite_face( i, B );
      Plane<T>    const plane    = make_plane(triangle);

      T const a0 = get_signed_distance(A0, plane);
      T const a1 = get_signed_distance(A1, plane);
      T const a2 = get_signed_distance(A2, plane);
      T const a3 = get_signed_distance(A3, plane);

      if ( a0 > 0 && a1 > 0 && a2 > 0 && a3 > 0 )
        return false;
    }

    // None of the face normals was a separating axis, now we test the edges

    // let edge0 be vertices 0 and 1
    // let edge1 be vertices 0 and 2
    // let edge2 be vertices 0 and 3
    // let edge3 be vertices 1 and 2
    // let edge4 be vertices 1 and 3
    // let edge5 be vertices 2 and 3

    // Ae0 v Be0
    const EigenVector3<T> ae0be0 = (A1-A0).cross( B1-B0);
    if (detail::SAT_line(A1, A2, A3, B1, B2, B3, ae0be0))
      return false;

    // Ae0 v Be1
    const EigenVector3<T> ae0be1 = (A1-A0).cross( B2-B0);
    if (detail::SAT_line(A1, A2, A3, B0, B1, B3, ae0be1))
      return false;

    // Ae0 v Be2
    const EigenVector3<T> ae0be2 = (A1-A0).cross( B3-B0);
    if (detail::SAT_line(A1, A2, A3, B0, B1, B2, ae0be2))
      return false;

    // Ae0 v Be3
    const EigenVector3<T> ae0be3 = (A1-A0).cross( B2-B1);
    if (detail::SAT_line(A1, A2, A3, B0, B1, B3, ae0be3))
      return false;

    // Ae0 v Be4
    const EigenVector3<T> ae0be4 = (A1-A0).cross( B3-B1);
    if (detail::SAT_line(A1, A2, A3, B0, B1, B2, ae0be4))
      return false;

    // Ae0 v Be5
    const EigenVector3<T> ae0be5 = (A1-A0).cross( B3-B2);
    if (detail::SAT_line(A1, A2, A3, B0, B1, B2, ae0be5))
      return false;

    // Ae1 v Be0
    const EigenVector3<T> ae1be0 = (A2-A0).cross( B1-B0);
    if (detail::SAT_line(A1, A2, A3, B1, B2, B3, ae1be0))
      return false;

    // Ae1 v Be1
    const EigenVector3<T> ae1be1 = (A2-A0).cross( B2-B0);
    if (detail::SAT_line(A1, A2, A3, B1, B2, B3, ae1be1))
      return false;

    // Ae1 v Be2
    const EigenVector3<T> ae1be2 = (A2-A0).cross( B3-B0);
    if (detail::SAT_line(A1, A2, A3, B0, B1, B2, ae1be2))
      return false;

    // Ae1 v Be3
    const EigenVector3<T> ae1be3 = (A2-A0).cross( B2-B1);
    if (detail::SAT_line(A1, A2, A3, B0, B1, B3, ae1be3))
      return false;

    // Ae1 v Be4
    const EigenVector3<T> ae1be4 = (A2-A0).cross( B3-B1);
    if (detail::SAT_line(A1, A2, A3, B0, B1, B2, ae1be4))
      return false;

    // Ae1 v Be5
    const EigenVector3<T> ae1be5 = (A2-A0).cross( B3-B2);
    if (detail::SAT_line(A1, A2, A3, B0, B1, B2, ae1be5))
      return false;

    // Ae2 v Be0
    const EigenVector3<T> ae2be0 = (A3-A0).cross( B1-B0);
    if (detail::SAT_line(A1, A2, A3, B1, B2, B3, ae2be0))
      return false;

    // Ae2 v Be1
    const EigenVector3<T> ae2be1 = (A3-A0).cross( B2-B0);
    if (detail::SAT_line(A1, A2, A3, B1, B2, B3, ae2be1))
      return false;

    // Ae2 v Be2
    const EigenVector3<T> ae2be2 = (A3-A0).cross( B3-B0);
    if (detail::SAT_line(A1, A2, A3, B1, B2, B3, ae2be2))
      return false;

    // Ae2 v Be3
    const EigenVector3<T> ae2be3 = (A3-A0).cross(B2-B1);
    if (detail::SAT_line(A1, A2, A3, B0, B1, B3, ae2be3))
      return false;

    // Ae2 v Be4
    const EigenVector3<T> ae2be4 = (A3-A0).cross( B3-B1);
    if (detail::SAT_line(A1, A2, A3, B0, B1, B2, ae2be4))
      return false;

    // Ae2 v Be5
    const EigenVector3<T> ae2be5 = (A3-A0).cross( B3-B2);
    if (detail::SAT_line(A1, A2, A3, B0, B1, B2, ae2be5))
      return false;

    // Ae3 v Be0
    const EigenVector3<T> ae3be0 = (A2-A1).cross(B1-B0);
    if (detail::SAT_line(A0, A1, A3, B1, B2, B3, ae3be0))
      return false;

    // Ae3 v Be1
    const EigenVector3<T> ae3be1 = (A2-A1).cross( B2-B0);
    if (detail::SAT_line(A0, A1, A3, B1, B2, B3, ae3be1))
      return false;

    // Ae3 v Be2
    const EigenVector3<T> ae3be2 = (A2-A1).cross( B3-B0);
    if (detail::SAT_line(A0, A1, A3, B1, B2, B3, ae3be2))
      return false;

    // Ae3 v Be3
    const EigenVector3<T> ae3be3 = (A2-A1).cross( B2-B1);
    if (detail::SAT_line(A0, A1, A3, B0, B1, B3, ae3be3))
      return false;

    // Ae3 v Be4
    const EigenVector3<T> ae3be4 = (A2-A1).cross( B3-B1);
    if (detail::SAT_line(A0, A1, A3, B0, B1, B2, ae3be4))
      return false;

    // Ae3 v Be5
    const EigenVector3<T> ae3be5 = (A2-A1).cross(B3-B2);
    if (detail::SAT_line(A0, A1, A3, B0, B1, B2, ae3be5))
      return false;

    // Ae4 v Be0
    const EigenVector3<T> ae4be0 = (A3-A1).cross( B1-B0);
    if (detail::SAT_line(A0, A1, A2, B1, B2, B3, ae4be0))
      return false;

    // Ae4 v Be1
    const EigenVector3<T> ae4be1 = (A3-A1).cross(B2-B0);
    if (detail::SAT_line(A0, A1, A2, B1, B2, B3, ae4be1))
      return false;

    // Ae4 v Be2
    const EigenVector3<T> ae4be2 = (A3-A1).cross(B3-B0);
    if (detail::SAT_line(A0, A1, A2, B1, B2, B3, ae4be2))
      return false;

    // Ae4 v Be3
    const EigenVector3<T> ae4be3 = (A3-A1).cross( B2-B1);
    if (detail::SAT_line(A0, A1, A2, B0, B1, B3, ae4be3))
      return false;

    // Ae4 v Be4
    const EigenVector3<T> ae4be4 = (A3-A1).cross(B3-B1);
    if (detail::SAT_line(A0, A1, A2, B0, B1, B2, ae4be4))
      return false;

    // Ae4 v Be5
    const EigenVector3<T> ae4be5 = (A3-A1).cross( B3-B2);
    if (detail::SAT_line(A0, A1, A2, B0, B1, B2, ae4be5))
      return false;

    // Ae5 v Be0
    const EigenVector3<T> ae5be0 = (A3-A2).cross( B1-B0);
    if (detail::SAT_line(A0, A1, A2, B1, B2, B3, ae5be0))
      return false;

    // Ae5 v Be1
    const EigenVector3<T> ae5be1 = (A3-A2).cross( B2-B0);
    if (detail::SAT_line(A0, A1, A2, B1, B2, B3, ae5be1))
      return false;

    // Ae5 v Be2
    const EigenVector3<T> ae5be2 = (A3-A2).cross( B3-B0);
    if (detail::SAT_line(A0, A1, A2, B1, B2, B3, ae5be2))
      return false;

    // Ae5 v Be3
    const EigenVector3<T> ae5be3 = (A3-A2).cross( B2-B1);
    if (detail::SAT_line(A0, A1, A2, B0, B1, B3, ae5be3))
      return false;

    // Ae5 v Be4
    const EigenVector3<T> ae5be4 = (A3-A2).cross( B3-B1);
    if (detail::SAT_line(A0, A1, A2, B0, B1, B2, ae5be4))
      return false;

    // Ae5 v Be5
    const EigenVector3<T> ae5be5 = (A3-A2).cross( B3-B2);
    if (detail::SAT_line(A0, A1, A2, B0, B1, B2, ae5be5))
      return false;

    return true;
  }

}// namespace geometry

// GEOMETRY_OVERLAP_TETRAHEDRON_TETRAHEDRON_H
#endif
