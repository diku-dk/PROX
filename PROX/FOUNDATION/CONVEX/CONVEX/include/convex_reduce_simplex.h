#ifndef CONVEX_REDUCE_SIMPLEX_H
#define CONVEX_REDUCE_SIMPLEX_H

#include <convex_reduce_edge.h>
#include <convex_reduce_triangle.h>
#include <convex_reduce_tetrahedron.h>

#include <cassert>

namespace convex
{
  /**
   * Reduce Simplex.
   * Compute the point, v, on the simplex that is closest to the origin and
   * Reduce simplex to lowest dimensional face on the boundary
   * containing the closest point.
   *
   * @param S         Upon invocation this argument holds the initial
   *                  simplex and upon return the argument holds the
   *                  resulting reduced simplex.
   *
   * @param a         Upon return this argument holds the corresponding
   *                  closest point on the set A.
   *
   * @param b         Upon return this argument holds the corresponding
   *                  closest point on the set B.
   *
   * @return          The closest point, v, on the simplex to the origin.
   */
  template<typename T>
inline EigenVector3<T> reduce_simplex( Simplex<T> & S, EigenVector3<T>& a, EigenVector3<T>& b )
  {

      const EigenVector3<T> p = {0, 0, 0};

    switch( dimension( S ) )
    {
      case 1:
      {
        // Nothing to do, a vertex can not be reduced!
        int    bit_A = 0;
        size_t idx_A = 0;

        get_used_indices( S.m_bitmask, idx_A, bit_A );

        S.m_bitmask = bit_A;        // 2011-11-12 Kenny: Why do we set the bitmask? It already has this value?

        S.m_w[idx_A] = 1;

        // 2011-11-12 Kenny: Are we certain that all other w's are zero? If not then we might get 'garbage' into the computation of the closest point? I think we should clear all w's and then set the idxA value to one!
      }
        break;
      case 2:
        reduce_edge(p, S);
        break;
      case 3:
        reduce_triangle(p, S);
        break;
      case 4:
        reduce_tetrahedron(p, S);
        break;
      default:
        assert(false || !"reduce_simplex(): can not reduce simplex of that size");
        break;
    };

    // Now compute the actual closest points based on the bary-centric coordinates.
    a = {0,0,0};
    b = {0,0,0};
    int used_bit = 1;
    for(size_t i=0u; i<4u; ++i)
    {
      if(S.m_bitmask & used_bit)
      {
        a += S.m_w[i]*S.m_a[i];
        b += S.m_w[i]*S.m_b[i];
      }
      used_bit <<= 1;
    }
    return (a-b);
  }

} // namespace convex

// CONVEX_REDUCE_SIMPLEX_H
#endif
