#ifndef CONVEX_REDUCE_TRIANGLE_H
#define CONVEX_REDUCE_TRIANGLE_H

#include <convex_outside_vertex_edge_voronoi_plane.h>
#include <convex_outside_edge_face_voronoi_plane.h>


#include <barycentric/geometry_barycentric.h>

#include <eigenhelperall.h>
#include <convex_simplex.h>
namespace convex
{

  /**
   * Reduce Triangle.
   * This function implements the case where a simplex
   * is a triangle. The function will compute the closest
   * point to p on the simplex and try to reduce the simplex
   * to the lowest dimensional face on the boundary of the
   * simplex containing the closest point.
   *
   * @param p   The test point.
   * @param S   Initially this argument holds the triangle simplex. Upon
   *            return the argument holds the reduced simplex.
   */
  template<typename T>
inline void reduce_triangle(const EigenVector3<T>& p_in, Simplex<T> & S)
  {

    int bit_A = 0;
    int bit_B = 0;
    int bit_C = 0;
    size_t idx_A = 0u;
    size_t idx_B = 0u;
    size_t idx_C = 0u;

    get_used_indices( S.m_bitmask, idx_A, bit_A, idx_B, bit_B, idx_C, bit_C );

    T scale = (S.m_v[idx_A]).norm() > (S.m_v[idx_B]).norm() ? (S.m_v[idx_A]).norm() : (S.m_v[idx_B]).norm();
    scale = scale > (S.m_v[idx_C]).norm() ? scale : (S.m_v[idx_C]).norm();

    const EigenVector3<T>& A = S.m_v[idx_A]/scale;
    const EigenVector3<T>& B = S.m_v[idx_B]/scale;
    const EigenVector3<T>& C = S.m_v[idx_C]/scale;
    const EigenVector3<T>& p = p_in/scale;

    bool const outside_AB  = outside_vertex_edge_voronoi_plane(p, A, B);
    bool const outside_AC  = outside_vertex_edge_voronoi_plane(p, A, C);
    bool const outside_BC  = outside_vertex_edge_voronoi_plane(p, B, C);
    bool const outside_BA  = outside_vertex_edge_voronoi_plane(p, B, A);
    bool const outside_CA  = outside_vertex_edge_voronoi_plane(p, C, A);
    bool const outside_CB  = outside_vertex_edge_voronoi_plane(p, C, B);
    bool const outside_ABC = outside_edge_face_voronoi_plane( p, A, B, C);
    bool const outside_BCA = outside_edge_face_voronoi_plane( p, B, C, A);
    bool const outside_CAB = outside_edge_face_voronoi_plane( p, C, A, B);

    // Test if we are in a vertex voronoi region
    if( outside_AB && outside_AC )
    {
      S.m_bitmask = bit_A;
      S.m_v[idx_B]= {0,0,0};
      S.m_v[idx_C]= {0,0,0};
      S.m_a[idx_B]= {0,0,0};
      S.m_a[idx_C]= {0,0,0};
      S.m_b[idx_B]= {0,0,0};
      S.m_b[idx_C]= {0,0,0};
      S.m_w[idx_A] = 1;
      S.m_w[idx_B] = 0;
      S.m_w[idx_C] = 0;
      return;
    }
    if( outside_BA && outside_BC )
    {
      S.m_bitmask = bit_B;
      S.m_v[idx_A]= {0,0,0};
      S.m_v[idx_C]= {0,0,0};
      S.m_a[idx_A]= {0,0,0};
      S.m_a[idx_C]= {0,0,0};
      S.m_b[idx_A]= {0,0,0};
      S.m_b[idx_C]= {0,0,0};
      S.m_w[idx_A] = 0;
      S.m_w[idx_B] = 1;
      S.m_w[idx_C] = 0;
      return;
    }
    if( outside_CA && outside_CB )
    {
      S.m_bitmask = bit_C;
      S.m_v[idx_A]= {0,0,0};
      S.m_v[idx_B]= {0,0,0};
      S.m_a[idx_A]= {0,0,0};
      S.m_a[idx_B]= {0,0,0};
      S.m_b[idx_A]= {0,0,0};
      S.m_b[idx_B]= {0,0,0};
      S.m_w[idx_A] = 0;
      S.m_w[idx_B] = 0;
      S.m_w[idx_C] = 1;
      return;
    }
    // Test voronoi regions of edges
    if( outside_ABC && !outside_AB && !outside_BA )
    {
      S.m_bitmask = bit_A | bit_B;
      S.m_v[idx_C]= {0,0,0};
      S.m_a[idx_C]= {0,0,0};
      S.m_b[idx_C]= {0,0,0};
      S.m_w[idx_C] = 0;
      geometry::barycentric(A,B,p,S.m_w[idx_A],S.m_w[idx_B]);
      return;
    }
    if( outside_BCA && !outside_BC && !outside_CB )
    {
      S.m_bitmask = bit_B | bit_C;
      S.m_v[idx_A]= {0,0,0};
      S.m_a[idx_A]= {0,0,0};
      S.m_b[idx_A]= {0,0,0};
      S.m_w[idx_A] = 0;
      geometry::barycentric(B,C,p,S.m_w[idx_B],S.m_w[idx_C]);
      return;
    }
    if( outside_CAB && !outside_AC && !outside_CA )
    {
      S.m_bitmask = bit_A | bit_C;
      S.m_v[idx_B]= {0,0,0};
      S.m_a[idx_B]= {0,0,0};
      S.m_b[idx_B]= {0,0,0};
      S.m_w[idx_B] = 0;
      geometry::barycentric(A,C,p,S.m_w[idx_A],S.m_w[idx_C]);
      return;
    }
    // Test if we are inside internal face voronnoi region
    if(!outside_ABC && !outside_BCA && !outside_CAB)
    {
      geometry::barycentric( A, B, C, p, S.m_w[idx_A], S.m_w[idx_B], S.m_w[idx_C]);
      return;
    }

  }

} // namespace convex

// CONVEX_REDUCE_TRIANGLE_H
#endif
