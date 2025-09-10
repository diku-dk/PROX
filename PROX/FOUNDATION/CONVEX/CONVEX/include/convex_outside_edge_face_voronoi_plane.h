#ifndef CONVEX_OUTSIDE_EDGE_FACE_VORONOI_PLANE_H
#define CONVEX_OUTSIDE_EDGE_FACE_VORONOI_PLANE_H


#include <cassert>
#include <eigenhelperall.h>

namespace convex
{

  /**
   * Test if Point is outside Edge Face Voronoi Plane.
   *
   * @param p     A point which we want to know whether it lies outside a face voronoi
   *              region. A face-voronoi region is given by a collection of edge-face voronoi
   *              planes and the face-plane.
   * @param A     The first vertex of the edge.
   * @param B     The second vertex of the edge.
   * @param C     The opposing vertex of the triangle containing the edge running from
   *              A to B. As such the C-vertex defines the plane with respect to which we want
   *              to test p.
   *
   * @return     If point is outside or on the plane then the return value is true otherwise it is false.
   */
  template< typename T >
  inline bool outside_edge_face_voronoi_plane(
                                       const EigenVector3<T>& p
                                       , const EigenVector3<T>& A
                                       , const EigenVector3<T>& B
                                       , const EigenVector3<T>& C
                                       )
  {


      const EigenVector3<T> m = ( A-C).cross(B-C );
    assert( dot( m, m ) > 0 || !"outside_edge_face_voronoi_plane(): Degenerate triangle encountered");

      const EigenVector3<T> n      = ( B-A).cross( m );
    T const sign_p = n.dot( p-B );
      T const sign_C = n.dot( C-B );

    assert( is_number( sign_p ) || !"outside_edge_face_voronoi_plane(): Not a Number encountered");
    assert( is_number( sign_C ) || !"outside_edge_face_voronoi_plane(): Not a Number encountered");

    return (sign_p*sign_C) <= 0;
  }

} // namespace convex

// CONVEX_OUTSIDE_EDGE_FACE_VORONOI_PLANE_H
#endif
