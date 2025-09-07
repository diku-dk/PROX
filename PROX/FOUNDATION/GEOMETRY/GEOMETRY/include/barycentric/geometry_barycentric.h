#ifndef GEOMETRY_BARYCENTRIC_H
#define GEOMETRY_BARYCENTRIC_H

#include <tiny_math_types.h>
#include <tiny_is_number.h>
#include <tiny_type_traits.h>
#include <tiny_vector.h>
#include <tiny_vector_functions.h>
#include <tiny_matrix_functions.h>

#include <cmath>
#include <cassert>

namespace geometry
{

  template<typename T>
  inline void barycentric(
      const EigenVector3<T>& x1
      , const EigenVector3<T>& x2
      , const EigenVector3<T> & p
      , T & w1
      , T& w2
      )
  {


      auto u  = x2-x1;
      T const uu = dot(u,u);

      assert( is_number(uu)   || !"barycentric(): NaN encountered");
      assert( uu > 0 || !"barycentric(): Degenerate edge encountered");

      // Project p onto edge running from x1 to x2.
      auto q  = (u.dot((p-x1))/ uu )*u + x1;

      auto a  = q - x2;
      auto b  = q - x1;

      T const aa =  dot(a,a) ;
      T const bb =  dot(b,b) ;

      assert( is_number(aa) || !"barycentric(): NaN encountered");
      assert( is_number(bb) || !"barycentric(): NaN encountered");

      w1 = ( u.dot(a) <= 0 ) ? std::sqrt( aa  / uu ) : - std::sqrt( aa  / uu );
      w2 = ( u.dot(b) >= 0 ) ? std::sqrt( bb  / uu ) : - std::sqrt( bb  / uu );

      assert( is_number(w1) || !"barycentric(): NaN encountered");
      assert( is_number(w2) || !"barycentric(): NaN encountered");
  }

  /*
   * Compute Barycentric Coordinates.
   * This method computes the barycentric coodinates for a point p of a triangle
   * given by the points x1,x2, and x3 (in counter clockwise order).
   *
   * The barycentric coordinates w1,w2, and w3 are defined such
   * that p' = w1*x1 + w2*x2 + w3*x3, is the point in plane of the
   * triangle closest to p.
   *
   * if 0 <= w1,w2,w3 <= 1 then the point lies inside or on the perimeter of the triangle.
   *
   * @warning  This method uses a geometric approach to compute the barycentric coordinates.
   *
   * @param x1    The first point of the triangle.
   * @param x2    The second point of the triangle.
   * @param x3    The third point of the triangle.
   * @param p     The point for which the barycentric coordinates should be computed.
   * @param w1    Upon return this parameter contains the value of the first barycentric coordinate.
   * @param w2    Upon return this parameter contains the value of the second barycentric coordinate.
   * @param w3    Upon return this parameter contains the value of the third barycentric coordinate.
   */
  template<typename T>
  inline void barycentric(
      const EigenVector3<T>& x1
      , const EigenVector3<T>& x2
      , const EigenVector3<T>& x3
      , const EigenVector3<T>& p
      , T& w1
      , T& w2
      , T& w3
      )
  {
      using std::sqrt;

      // Before computing the barycentric coordinates, we 'normalize' the triangle to ensure it is not a sliver
      auto b1 = x2 - x1;
      auto b2 = x3 - x1;
      auto b3 = ((b1).normalized()).cross((b2).normalized());

      const EigenMatrix3<T> basis{{b1(0), b2(0), b3(0)},
                                  {b1(1), b2(1), b3(1)},
                                  {b1(2), b2(2), b3(2)}};

/*      M const basis = M::make(  b1(0), b2(0), b3(0)
                              , b1(1), b2(1), b3(1)
                              , b1(2), b2(2), b3(2));*/

      // q is the point p transformed to the isoparametric system
      const EigenVector3<T> q = basis.inverse()*(p - x1);

      w1 = 1 - q[0] - q[1];
      w2 = q[0];
      w3 = q[1];

      assert( is_number(w1) || !"barycentric(): NaN encountered");
      assert( is_number(w2) || !"barycentric(): NaN encountered");
      assert( is_number(w3) || !"barycentric(): NaN encountered");
  }

  /*
   * Compute Barycentric Coordinates.
   * This method computes the barycentric coodinates for a point p of a tetrahedron
   * given by the points x1,x2,x3, x4 (in right-hand order).
   *
   * @warning  This method uses a geometric approach to compute the barycentric coordinates.
   *
   * @param x1    The first point of the triangle.
   * @param x2    The second point of the triangle.
   * @param x3    The third point of the triangle.
   * @param x4    The fourth point of the triangle.
   * @param p     The point for which the barycentric coordinates should be computed.
   * @param w1    Upon return this parameter contains the value of the first barycentric coordinate.
   * @param w2    Upon return this parameter contains the value of the second barycentric coordinate.
   * @param w3    Upon return this parameter contains the value of the third barycentric coordinate.
   * @param w4    Upon return this parameter contains the value of the fourth barycentric coordinate.
   */
  template<typename T>
  inline void barycentric(
      const EigenVector3<T>& x1
      , const EigenVector3<T>& x2
      , const EigenVector3<T>& x3
      , const EigenVector3<T>& x4
      , const EigenVector3<T>& p
      , T& w1
      , T& w2
      , T& w3
      , T& w4
      )
  {
      using std::fabs;

      const EigenVector3<T> b1 = x2 - x1;
      const EigenVector3<T> b2 = x3 - x1;
      const EigenVector3<T> b3 = x4 - x1;

      // check if basis forms a right hand coordsys, if not exchange two base vectors
      if(dot(b3,(b1.cross(b2))) > 0)
      {
          const EigenMatrix3<T> basis{{b1(0), b2(0), b3(0)},
              {b1(1), b2(1), b3(1)},
              {b1(2), b2(2), b3(2)}};
          // q is the point p transformed to the isoparametric system
          const EigenVector3<T> q = (basis).inverse()*(p - x1);

          w1 = 1 - q[0] - q[1] - q[2];
          w2 = q[0];
          w3 = q[1];
          w4 = q[2];
      }
      else
      {
          const EigenMatrix3<T> basis{{b1(0), b2(0), b3(0)},
              {b1(1), b2(1), b3(1)},
              {b1(2), b2(2), b3(2)}};
          // q is the point p transformed to the isoparametric system
          const EigenVector3<T> q = (basis).inverse()*(p - x1);

          w1 = 1 - q[0] - q[1] - q[2];
          w2 = q[0];
          w3 = q[2];// we alter the sequence of barycentric coords to match altered basis
          w4 = q[1];
      }

      assert( is_number(w1) || !"barycentric(): NaN encountered");
      assert( is_number(w2) || !"barycentric(): NaN encountered");
      assert( is_number(w3) || !"barycentric(): NaN encountered");
      assert( is_number(w4) || !"barycentric(): NaN encountered");
  }

}//namespace geometry

// GEOMETRY_BARYCENTRIC_H
#endif
