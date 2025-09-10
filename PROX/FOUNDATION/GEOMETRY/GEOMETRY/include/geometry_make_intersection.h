#ifndef GEOMETRY_MAKE_INTERSECTION_H
#define GEOMETRY_MAKE_INTERSECTION_H

#include <types/geometry_line.h>
#include <types/geometry_plane.h>

#include <eigenhelperall.h>

#include <cassert>
#include <cmath>  // needed for std::fabs

namespace geometry
{

template <typename T>
inline EigenVector3<T> make_intersection(Line<T> const& L, Plane<T> const& P)
{
    using std::fabs;
    const EigenVector3<T>& o       = L.point();
    const EigenVector3<T>& d       = L.direction();
    const EigenVector3<T>& n       = P.normal();
    T const & w       = P.offset();

    //
    // A line is given as all points
    //
    //     p(t) = o + t*d for all t
    //
    // For all points on the line that are on the plane we have
    //
    //     n*p(t) - w = 0
    //
    // Substitution yiedls
    //
    //     n*o + t*n*d - w = 0
    //
    //  we find t-value as
    //
    //     t =  ( w - n*o ) / n*d
    //
    //

    T const   n_dot_d = dot( n, d );

    assert( is_number(n_dot_d)         || !"make_intersection(): NaN encountered");
    assert( is_finite(n_dot_d)         || !"make_intersection(): Inf encountered");
    assert(fabs( n_dot_d) > 0 || !"make_intersection(): line was parallel withe plane");


    T const   n_dot_o = dot( n, o );

    assert( is_number(n_dot_o)         || !"make_intersection(): NaN encountered");
    assert( is_finite(n_dot_o)         || !"make_intersection(): Inf encountered");

    T const   t       = (w - n_dot_o) / n_dot_d;

    assert( is_number(t)               || !"make_intersection(): NaN encountered");
    assert( is_finite(t)               || !"make_intersection(): Inf encountered");
    const  EigenVector3<T> p       = o + t*d;

    assert( is_number(p(0))               || !"make_intersection(): NaN encountered");
    assert( is_finite(p(0))               || !"make_intersection(): Inf encountered");
    assert( is_number(p(1))               || !"make_intersection(): NaN encountered");
    assert( is_finite(p(1))               || !"make_intersection(): Inf encountered");
    assert( is_number(p(2))               || !"make_intersection(): NaN encountered");
    assert( is_finite(p(2))               || !"make_intersection(): Inf encountered");

    return p;
  }

  template<typename T>
  inline Line<T> make_intersection( Plane<T> const & A,  Plane<T> const & B)
  {
    const EigenVector3<T> D = cross( A.normal(), B.normal() );

    assert(dot(D, D) > 0 || !"make_intersection(): error planes are coplanar");
    assert( is_number(D(0))               || !"make_intersection(): NaN encountered");
    assert( is_finite(D(0))               || !"make_intersection(): Inf encountered");
    assert( is_number(D(1))               || !"make_intersection(): NaN encountered");
    assert( is_finite(D(1))               || !"make_intersection(): Inf encountered");
    assert( is_number(D(2))               || !"make_intersection(): NaN encountered");
    assert( is_finite(D(2))               || !"make_intersection(): Inf encountered");

    const EigenVector3<T> d = unit( D );

    assert( is_number(d(0))               || !"make_intersection(): NaN encountered");
    assert( is_finite(d(0))               || !"make_intersection(): Inf encountered");
    assert( is_number(d(1))               || !"make_intersection(): NaN encountered");
    assert( is_finite(d(1))               || !"make_intersection(): Inf encountered");
    assert( is_number(d(2))               || !"make_intersection(): NaN encountered");
    assert( is_finite(d(2))               || !"make_intersection(): Inf encountered");

    // Get a point on plane A
    const EigenVector3<T> pA = A.normal()*A.offset();

    assert( is_number(pA(0))               || !"make_intersection(): NaN encountered");
    assert( is_finite(pA(0))               || !"make_intersection(): Inf encountered");
    assert( is_number(pA(1))               || !"make_intersection(): NaN encountered");
    assert( is_finite(pA(1))               || !"make_intersection(): Inf encountered");
    assert( is_number(pA(2))               || !"make_intersection(): NaN encountered");
    assert( is_finite(pA(2))               || !"make_intersection(): Inf encountered");

    // Direction on plane A towards plane B
    const EigenVector3<T> r = unit(cross( A.normal(), d ));

    assert( is_number(r(0))               || !"make_intersection(): NaN encountered");
    assert( is_finite(r(0))               || !"make_intersection(): Inf encountered");
    assert( is_number(r(1))               || !"make_intersection(): NaN encountered");
    assert( is_finite(r(1))               || !"make_intersection(): Inf encountered");
    assert( is_number(r(2))               || !"make_intersection(): NaN encountered");
    assert( is_finite(r(2))               || !"make_intersection(): Inf encountered");

    // Now we find the intersection point between the line on plane A with origin pA and direction r and plane B
    //
    // Since line (pA, r) is on plane A then the intersection point with B will be a point both on plane A and B.
    //
    Line<T> const L(pA, r);

    const EigenVector3<T> o = make_intersection(L, B);

    assert( is_number(o(0))               || !"make_intersection(): NaN encountered");
    assert( is_finite(o(0))               || !"make_intersection(): Inf encountered");
    assert( is_number(o(1))               || !"make_intersection(): NaN encountered");
    assert( is_finite(o(1))               || !"make_intersection(): Inf encountered");
    assert( is_number(o(2))               || !"make_intersection(): NaN encountered");
    assert( is_finite(o(2))               || !"make_intersection(): Inf encountered");

    return make_line(o, d, FROM_DIRECTION() );
  }

  template<typename T>
  inline EigenVector3<T> make_intersection(Plane<T> const & A, Plane<T> const & B, Plane<T> const & C)
  {
    Line<T> const L = make_intersection(A,B);
    return make_intersection(L,C);
  }

}// namespace geometry

// GEOMETRY_MAKE_INTERSECTION_H
#endif
