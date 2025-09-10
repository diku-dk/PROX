#ifndef CONVEX_OUTSIDE_TRIANGLE_H
#define CONVEX_OUTSIDE_TRIANGLE_H

#include <cassert>
#include <eigenhelperall.h>

namespace convex
{

  /**
   * Test if point is outside a Triangle.
   *
   * @param p     A point which we want to know whether it lies on the other
   *              side of the triangle that the point q.
   * @param A     The first vertex of the triangle.
   * @param B     The second vertex of the triangle.
   * @param C     The third vertex of the triangle.
   * @param q     A point that is known to lie on the back-side of the triangle
   *
   *@return       If p is outside or on the face plane then the return value is true otherwise it is false.
   */
  template< typename T>
  inline bool outside_triangle(
                        const EigenVector3<T>& p
                        , const EigenVector3<T>& A
                        , const EigenVector3<T>& B
                        , const EigenVector3<T>& C
                        , const EigenVector3<T>& q
                        )
  {
      const EigenVector3<T> n = ( A-B).cross( C-B );

    assert( dot( n, n ) > 0 || !"outside_triangle(): Degenerate triangle encountered");

      T const sign_p = ( n).dot( p-B );
    T const sign_q = ( n).dot( q-B );

    assert( is_number( sign_p ) || !"outside_triangle(): Not a Number encountered");
    assert( is_number( sign_q ) || !"outside_triangle(): Not a Number encountered");

    assert( sign_q < 0 || sign_q > 0 || !"outside_triangle(): q was in plane, can  not be used to determine sign");

    return (sign_p*sign_q) <= 0;
  }

} // namespace convex

// CONVEX_OUTSIDE_TRIANGLE_H
#endif
