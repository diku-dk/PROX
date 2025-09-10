#ifndef GEOMETRY_BARYCENTRIC_DERIVATIVES_H
#define GEOMETRY_BARYCENTRIC_DERIVATIVES_H

#include <eigenhelperall.h>
#include <cassert>

namespace geometry
{

  template<typename T>
  inline void barycentric_derivatives(
                                      const EigenVector3<T>& x1
                                      , const EigenVector3<T>& x2
                                      , const EigenVector3<T>& x3
                                      , const EigenVector3<T>& x4
                                      , EigenVector3<T>& nabla_w1
                                      , EigenVector3<T>& nabla_w2
                                      , EigenVector3<T>& nabla_w3
                                      , EigenVector3<T>& nabla_w4
                                      )
  {
    using std::fabs;

    T const vol6 = dot((x4 - x1).eval(), (cross(x2 - x1, x3 - x1)));

    assert( vol6 > 0 || !"compute_barycentric_derivatives(): tetrahedron was left oriented");

    nabla_w1 = cross( x4 - x2, x3 - x2) / vol6;
    nabla_w2 = cross( x3 - x1, x4 - x1) / vol6;
    nabla_w3 = cross( x4 - x1, x2 - x1) / vol6;
    nabla_w4 = cross( x2 - x1, x3 - x1) / vol6;
  }

}//namespace geometry

// GEOMETRY_BARYCENTRIC_DERIVATIVES_H
#endif
