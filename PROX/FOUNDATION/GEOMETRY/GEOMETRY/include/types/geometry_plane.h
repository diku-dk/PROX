#ifndef GEOMETRY_PLANE_H
#define GEOMETRY_PLANE_H

#include <types/geometry_triangle.h>

#include <tiny_precision.h>
#include <tiny_math_types.h>

#include <cmath>


namespace geometry
{

  /**
   * Set of all points p such that
   *
   *  n^T (p-p0) = 0
   *
   * where p0 is a point in the plane and n is the outward plane normal.
   *
   * We define w as
   *
   *  w = n*p0
   *
   * Given that n is the unit normal then w is the distance form origin to the plane.
   */
  template<typename T>
  class Plane
  {
  protected:
      EigenVector3<T> m_normal;
      T m_offset;

  public:

    void set_normal(const EigenVector3<T>& normal)
    {
      this->m_normal = unit(normal);
    }

    T       & offset()       { return this->m_offset; }
    const EigenVector3<T>& normal() const { return this->m_normal; }
    T const & offset() const { return this->m_offset; }

    void set_n(const EigenVector3<T>& n)
    {
      this->m_normal = unit(n);
    }

    T       & w()       { return this->m_offset; }
    const EigenVector3<T>& n() const { return this->m_normal; }
    T const & w() const { return this->m_offset; }

  public:

    Plane()
          : m_normal( EigenVector3<T>(0,0,1) )
    , m_offset( 0 )
    {}

    ~Plane(){}

    Plane(const EigenVector3<T>& normal, T const & offset)
    : m_normal( unit(normal) )
    , m_offset(offset)
    {
      assert(fabs(1 - norm(normal)) < tiny::working_precision<T>() || !"Plane(): Must be unit normal");
    }

    Plane(Plane const & plane)
    {
      *this = plane;
    }

    Plane& operator=(Plane const & plane)
    {
      if( this != &plane)
      {
        this->m_normal = plane.m_normal;
        this->m_offset = plane.m_offset;
      }
      return *this;
    }


  };

  template<typename T>
  inline Plane<T> make_plane(const EigenVector3<T>& normal, const T& offset)
  {

    assert(fabs(1 - norm(normal)) < tiny::working_precision<T>() || !"make_plane(): Must be unit normal");

    return Plane<T>(normal,offset);
  }

  template<typename T>
  inline Plane<T> make_plane(const EigenVector3<T>& p0, const EigenVector3<T>& p1, const EigenVector3<T>& p2)
  {

      const EigenVector3<T> normal = ( (p1-p0).cross( p2-p0) ).normalized();
    T const offset = dot(p0,normal);

    return Plane<T>(normal,offset);
  }

  /**
   * Tries to counter numerical precision errors as much as possible.
   * This is 3 times more expensive to do that just making a plane.
   */
  template<typename T>
  inline Plane<T> make_precise_plane(const EigenVector3<T>& p0, const EigenVector3<T>& p1, const EigenVector3<T>& p2)
  {

    const EigenVector3<T> m0 = (p1-p0).cross (p2-p0);
    const EigenVector3<T> m1 = (p2-p1).cross (p0-p1);
    const EigenVector3<T> m2 = (p0-p2).cross( p1-p2);

    T const l0 = dot(m0,m0);
    T const l1 = dot(m1,m1);
    T const l2 = dot(m2,m2);

    if (l0 >= l1 && l0 >= l2)
    {
      const EigenVector3<T> normal = m0 / sqrt(l0);
      T const offset = inner_prod(p0,normal);

      return Plane<T>(normal,offset);
    }
    if (l1 >= l0 && l1 >= l2)
    {
      const EigenVector3<T> normal = m1 / sqrt(l1);
      T const offset = inner_prod(p1,normal);

      return Plane<T>(normal,offset);
    }
    if (l2 >= l1 && l2 >= l0)
    {
      const EigenVector3<T> normal = m2 / sqrt(l2);
      T const offset = inner_prod(p2,normal);

      return Plane<T>(normal,offset);
    }

    return Plane<T>();
  }

  template<typename T>
  inline Plane<T> make_plane(const EigenVector3<T>& normal, const EigenVector3<T>& p)
  {

    T const offset = dot(p ,normal);

    return Plane<T>(normal,offset);
  }

  template<typename T>
  inline Plane<T> make_plane(Triangle<T> const & t)
  {
    return make_plane<T>(t.point(0),t.point(1),t.point(2));
  }

  template<typename T>
  inline T get_signed_distance(const EigenVector3<T>& p, Plane<T>  const & P)
  {
    return dot(P.n(),p) - P.w();
  }

  template<typename T>
  inline T get_distance(const EigenVector3<T>& p, Plane<T>  const & P)
  {
    using std::fabs;

    return fabs( get_signed_distance(p,P));
  }



}// namespace geometry

// GEOMETRY_PLANE_H
#endif
