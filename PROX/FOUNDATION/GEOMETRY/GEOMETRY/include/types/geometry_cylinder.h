#ifndef GEOMETRY_CYLINDER_H
#define GEOMETRY_CYLINDER_H

#include <eigenhelperall.h>

namespace geometry
{

template <typename T> class CylinderEigen
{
protected:
    T m_radius;
    T m_height;
    EigenVector3<T> m_axis;
    EigenVector3<T> m_center;

public:
    T const& radius() const { return this->m_radius; }
    T const& height() const { return this->m_height; }
    T half_height() const { return this->m_height * 0.5f; }
    const EigenVector3<T>& axis() const { return this->m_axis; }
    const EigenVector3<T>& center() const { return this->m_center; }

    void set_radius(T const& radius)
    {
        assert(radius >= 0
               || !"CylinderEigen::set_radius(): Value was negative");

        this->m_radius = radius;
    }

    void set_height(T const& height)
    {
        assert(height >= 0
               || !"CylinderEigen::set_height(): Value was negative");

        this->m_height = height;
    }

    void set_axis(const EigenVector3<T>& axis)
    {
        assert(norm(axis) > 0
               || !"CylinderEigen::set_axis(): Axis was zero-vector");

        this->m_axis = unit(axis);
    }

    EigenVector3<T>& center() { return this->m_center; }

public:
    CylinderEigen()
        : m_radius(1)
        , m_height(1)
        , m_axis(EigenVector3<T>(0, 0, 1))
        , m_center(EigenVector3<T>(0, 0, 0))
    {
    }

    CylinderEigen(T const& radius, T const& height, const EigenVector3<T>& axis,
                  const EigenVector3<T>& center)
        : m_radius(radius)
        , m_height(height)
        , m_axis(unit(axis))
        , m_center(center)
    {
        assert(radius >= 0 || !"CylinderEigen(): radius was negative");
        assert(height >= 0 || !"CylinderEigen(): height was negative");
        assert(norm(axis) > 0 || !"CylinderEigen(): axis was zero-vector");
    }

    CylinderEigen(CylinderEigen const& CylinderEigen) { *this = CylinderEigen; }

    CylinderEigen& operator=(CylinderEigen const& CylinderEigen)
    {
        if (this != &CylinderEigen)
        {
            this->m_radius = CylinderEigen.m_radius;
            this->m_height = CylinderEigen.m_height;
            this->m_axis = CylinderEigen.m_axis;
            this->m_center = CylinderEigen.m_center;
        }
        return *this;
    }
  };

  template<typename T>
  CylinderEigen<T> make_cylinder(
      const T& radius
      , const T& height
      , const EigenVector3<T>& axis
      , const EigenVector3<T>& center
      )
  {
      return CylinderEigen<T>(radius, height, axis, center);
  }

  } //namespace geometry

// GEOMETRY_CYLINDER_H
#endif
