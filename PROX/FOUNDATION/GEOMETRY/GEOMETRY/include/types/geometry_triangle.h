#ifndef GEOMETRY_TRIANGLE_H
#define GEOMETRY_TRIANGLE_H

#include <eigenhelperall.h>

namespace geometry
{

  template<typename T>
  class Triangle
  {
  public:

  protected:
      EigenVector3<T> m_point[3];

  public:

    const EigenVector3<T>& point(unsigned int const & idx ) const { return m_point[idx]; }
    EigenVector3<T>& point(unsigned int const & idx )       { return m_point[idx]; }
    const EigenVector3<T>& p(unsigned int const & idx )     const { return m_point[idx]; }
    EigenVector3<T>& p(unsigned int const & idx )           { return m_point[idx]; }

  public:

    Triangle()
    {
        m_point[0] = {0,0,0};
        m_point[1] = {1,0,0};
        m_point[2] = {0,1,0};
    }

    Triangle(
                const EigenVector3<T>& p0
                , const EigenVector3<T>& p1
                , const EigenVector3<T>& p2
                )
    {
      m_point[0] = p0;
      m_point[1] = p1;
      m_point[2] = p2;
    }

    Triangle( Triangle const & tri)
    {
      *this = tri;
    }

    Triangle & operator=(Triangle const & tri)
    {
      if (this != &tri)
      {
        this->m_point[0] = tri.m_point[0];
        this->m_point[1] = tri.m_point[1];
        this->m_point[2] = tri.m_point[2];
      }
      return *this;
    }

  };

  template<typename T>
  inline Triangle<T> make_triangle()
  {
    return Triangle<T>();
  }

  template<typename T>
  inline Triangle<T> make_triangle(
                  const EigenVector3<T>& p0
                  , const EigenVector3<T>& p1
                  , const EigenVector3<T>& p2
                  )
  {
    return Triangle<T>(p0,p1,p2);
  }

  template<typename T>
  inline T get_area(Triangle<T> const & tri)
  {

    const EigenVector3<T>& x1 = tri.point(0);
    const EigenVector3<T>& x2 = tri.point(1);
    const EigenVector3<T>& x3 = tri.point(2);

    T const area2 = ( (x2 - x1).cross(x3 - x1) ).norm();

    assert( is_number(area2) || !"get_area(): Nan encountered");
    assert( is_finite(area2) || !"get_area(): Inf encountered");

    return area2 / 2;
  }

}// namespace geometry

//GEOMETRY_TRIANGLE_H
#endif
