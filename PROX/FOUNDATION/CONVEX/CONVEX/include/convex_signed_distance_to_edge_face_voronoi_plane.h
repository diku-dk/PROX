#ifndef CONVEX_SIGNED_DISTANCE_TO_EDGE_FACE_VORONOI_PLANE_H
#define CONVEX_SIGNED_DISTANCE_TO_EDGE_FACE_VORONOI_PLANE_H

#include <tiny_vector_functions.h>
#include <tiny_is_number.h>

#include <cmath>
#include <cassert>
#include <tiny_math_types.h>

namespace convex
{

  /**
   * Signed Distance between a point and a Edge Face Voronoi Plane.
   *
   * @param p A point which we want to know whether it lies outside a face voronoi
   * region. A face-voronoi region is given by a collection of edge-face voronoi
   * planes and the face-plane.
   * @param A    The first vertex of the edge.
   * @param B    The second vertex of the edge.
   * @param C The opposing vertex of the triangle containing the edge running from
   * A to B. As such the C-vertex defines the plane with respect to which we want
   * to test p.
   *
   * @return     The signed distance of p to the edge-face voronoi plane.
   */
  template< typename T>
  inline T signed_distance_to_edge_face_voronoi_plane(
                                                                   const EigenVector3<T>& p
                                                                   , const EigenVector3<T> & A
                                                                   , const EigenVector3<T> & B
                                                                   , const EigenVector3<T> & C
                                                                   )
  {
    using std::fabs;

      EigenVector3<T>  m      = ( A-C).cross(B-C );

    assert( tiny::inner_prod( m, m ) > 0 || !"signed_distance_to_edge_face_voronoi_plane(): Degenerate triangle encountered");

      EigenVector3<T>  l      = ( B-A).cross( m );
    EigenVector3<T>  n      = ( l ).norm();

      T sign_p = ( n).dot( p-B );
    T sign_C = ( n).dot( C-B );
    T abs_p  = fabs( sign_p );

    assert( is_number( sign_p ) || !"signed_distance_to_edge_face_voronoi_plane(): Not a Number encountered");
    assert( is_number( sign_C ) || !"signed_distance_to_edge_face_voronoi_plane(): Not a Number encountered");
    assert( is_number( abs_p )  || !"signed_distance_to_edge_face_voronoi_plane(): Not a Number encountered");

    bool in_front = ( (sign_p*sign_C) <= 0 );

    return in_front ? abs_p : -abs_p;
  }

} // namespace convex

// CONVEX_SIGNED_DISTANCE_TO_EDGE_FACE_VORONOI_PLANE_H
#endif
