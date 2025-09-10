#ifndef CONVEX_OUTSIDE_VERTEX_EDGE_VORONOI_PLANE_H
#define CONVEX_OUTSIDE_VERTEX_EDGE_VORONOI_PLANE_H

#include <eigenhelperall.h>
#include <cassert>

namespace convex
{

  /**
   * Test if point is outside a vertex edge voronoi plane.
   * The vertex edge voronoi plane is defined such that the plane normal is given
   * as the unit vector of the vector \fA-B\f and the point \fA\f is defined to
   * lie in the plane, and the point \fB\f behind the plane.
   *
   * @param p    The point that should be tested.
   * @param A    The first vertex of the edge.
   * @param B    The second vertex of the edge.
   *
   * @return     If point is outside or on plane then return value is true otherwise it is false.
   */
  template< typename T>
  inline bool outside_vertex_edge_voronoi_plane(
                                         const EigenVector3<T>& p
                                         , const EigenVector3<T>& A
                                         , const EigenVector3<T>& B
                                         )
  {

    EigenVector3<T> n = (A-B);
    assert( dot( n, n ) > 0 || !"outside_vertex_edge_voronoi_plane(): Degenerate edge encountered");

    T const sign_p = n.dot( (p-A) );
    assert( is_number( sign_p ) || !"outside_vertex_edge_voronoi_plane(): Not a Number encountered");

    return sign_p >= 0;
  }

} // namespace convex

// CONVEX_OUTSIDE_VERTEX_EDGE_VORONOI_PLANE_H
#endif
