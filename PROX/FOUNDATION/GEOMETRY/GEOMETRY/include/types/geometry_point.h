#ifndef GEOMETRY_POINT_H
#define GEOMETRY_POINT_H

#include <types/geometry_support_mapping.h>

namespace geometry
{

  template<typename V>
  class Point
    : public geometry::SupportMapping<typename V::real_type>
  {
  public:
      using T = typename V::real_type;
      using VT = typename V::value_traits;

  protected:
      V m_coordinates;

  public:

    V const & coord() const   { return this->m_coordinates; }
    V       & coord()         { return this->m_coordinates; }

  public:

    Point()
    : m_coordinates( V::zero() )
    {}

    Point(V const &  coordinates )
    : m_coordinates( coordinates )
    {}

    Point(Point const & p)
    {
      *this = p;
    }

    Point & operator=(Point const & p)
    {
      if( this != &p)
      {
        this->m_coordinates = p.m_coordinates;
      }
      return *this;
    }

  public:

    V get_support_point(V const & v) const
    {
      return this->m_coordinates;
    }
    EigenVector3<T> get_support_point(EigenVector3<T> dir) const override
    {
        return toEigen(m_coordinates);
    }


    T get_scale() const
    {
      return 0;
    }
  };

} //namespace geometry

// GEOMETRY_POINT_H
#endif
