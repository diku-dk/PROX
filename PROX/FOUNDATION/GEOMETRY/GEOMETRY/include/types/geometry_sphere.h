#ifndef GEOMETRY_SPHERE_H
#define GEOMETRY_SPHERE_H

#include <types/geometry_support_mapping.h>

#include <tiny_is_number.h>
#include <tiny_is_finite.h>

#include <cmath> // needed for std::sqrt
#include <cassert>

namespace geometry
{

  template<typename V>
  class Sphere
  : public geometry::SupportMapping<V>

  {
  public:
      using T = typename V::real_type;
      using VT = typename V::value_traits;

  protected:
      V m_center;
      T m_radius;

  public:

    V const & center() const { return this->m_center; }
    T const & radius() const { return this->m_radius; }
    V       & center()       { return this->m_center; }
    T       & radius()       { return this->m_radius; }

  public:

    Sphere()
      : m_center( V::zero() )
      , m_radius( 1 )
    {
    }

    Sphere( V const & center, T const &  radius )
      : m_center( center )
      , m_radius( radius )
    {
    }

    Sphere(Sphere const & sph)
    {
      *this = sph;
    }

    Sphere & operator=(Sphere const & sph)
    {
      if( this != &sph)
      {
        this->m_center = sph.m_center;
        this->m_radius = sph.m_radius;
      }
      return *this;
    }


  public:
    EigenVector3<T> get_support_point(EigenVector3<T> dir) const
    {
        auto dirLen = dir.norm();
        assert(is_number(dirLen) && is_finite(dirLen));
        if (dirLen > 0)
        {
            return (m_radius / dirLen) * dir + toEigen(m_center);
        }
        return toEigen(m_center) + EigenVector3<T>{m_radius, 0, 0};
    }

    T get_scale() const
    {
      assert( is_number(this->m_radius)   || !"get_scale(): NAN encountered");
      assert( is_finite(this->m_radius)   || !"get_scale(): INF encountered");
      assert( this->m_radius > 0 || !"get_scale(): radius was non-positive");

      return 2 *  this->m_radius ;
    }

  };

  template<typename V>
  inline Sphere<V> make_sphere(V const & center, typename V::real_type const &  radius)
  {
    return Sphere<V>(center,radius);
  }

}// namespace geometry

//GEOMETRY_SPHERE_H
#endif
