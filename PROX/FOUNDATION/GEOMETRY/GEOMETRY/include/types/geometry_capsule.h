#ifndef GEOMETRY_CAPSULE_H
#define GEOMETRY_CAPSULE_H

#include <eigenhelperall.h>

namespace geometry
{

template <typename T> class CapsuleEigen
{
protected:
    T m_radius;
    EigenVector3<T> m_point0;
    EigenVector3<T> m_point1;

public:
    T const& radius() const { return this->m_radius; }
    const EigenVector3<T>& point0() const { return this->m_point0; }
    const EigenVector3<T>& point1() const { return this->m_point1; }
    EigenVector3<T>& point0() { return this->m_point0; }
    EigenVector3<T>& point1() { return this->m_point1; }

    void set_radius(T const& radius)
    {
        assert(radius >= 0 || !"Capsule::set_radius(): Value was negative");

        this->m_radius = radius;
    }

public:
    CapsuleEigen()
        : m_radius(1)
        , m_point0(EigenVector3<T>(0, 0, 0))
        , m_point1(EigenVector3<T>(0, 0, 0))
    {
    }

    CapsuleEigen(T const& radius, const EigenVector3<T>& point0,
                 const EigenVector3<T>& point1)
        : m_radius(radius)
        , m_point0(point0)
        , m_point1(point1)
    {
        assert(radius >= 0 || !"Capsule(): radius was negative");
    }

    CapsuleEigen(CapsuleEigen const& capsule) { *this = capsule; }

    CapsuleEigen& operator=(CapsuleEigen const& capsule)
    {
        if (this != &capsule)
        {
            this->m_radius = capsule.m_radius;
            this->m_point0 = capsule.m_point0;
            this->m_point1 = capsule.m_point1;
        }
        return *this;
    }
  };

  template<typename T>
  CapsuleEigen<T> make_capsule(
      const T& radius
      , const EigenVector3<T>& point0
      , const EigenVector3<T>& point1
      )
  {
      return CapsuleEigen<T>(radius, point0, point1);
  }

  }

// GEOMETRY_CAPSULE_H
#endif
