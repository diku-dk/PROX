#ifndef CONVEX_SIGNED_DISTANCE_TO_VERTEX_EDGE_VORONOI_PLANE_H
#define CONVEX_SIGNED_DISTANCE_TO_VERTEX_EDGE_VORONOI_PLANE_H

#include <cassert>
#include <eigenhelperall.h>

namespace convex
{
  /**
   * Signed Distance between a point and a vertex edge voronoi plane.
   * The vertex edge voronoi plane is defined such that the plane normal is given
   * as the unit vector of the vector \fA-B\f and the point \fA\f is defined to
   * lie in the plane, and the point \fB\f behind the plane.
   *
   * @param p    The point that should be tested.
   * @param A    The first vertex of the edge.
   * @param B    The second vertex of the edge.
   *
   * @return     The signed distance of the point p.
   */
  template< typename T>
  inline T signed_distance_to_vertex_edge_voronoi_plane(const EigenVector3<T> & p, const EigenVector3<T> & A,
                                                        const EigenVector3<T> & B)
  {


    const EigenVector3<T> m = (A-B);

    assert( dot( m, m ) > 0 || !"signed_distance_to_vertex_edge_voronoi_plane(): Degenerate edge encountered");

    const EigenVector3<T> n = (m).normalized();

    T sign_p = (n).dot((p-A) );

    assert( is_number( sign_p ) || !"signed_distance_to_vertex_edge_voronoi_plane(): Not a Number encountered");

    return sign_p;
  }

} // namespace convex

// CONVEX_SIGNED_DISTANCE_TO_VERTEX_EDGE_VORONOI_PLANE_H
#endif
