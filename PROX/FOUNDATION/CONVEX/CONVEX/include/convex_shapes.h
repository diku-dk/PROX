#ifndef CONVEX_SHAPES_H
#define CONVEX_SHAPES_H

#include <types/geometry_support_mapping.h>

#include <vector>

namespace convex
{

  /**
   * The cylinder axis is by default equal to the z-axis.
   */
  template<typename T>
  class Cylinder
  : public geometry::SupportMapping<T>
  {

  protected:

    T    m_half_height;     ///< The half height of the cylinder. Default value is one.
    T    m_radius;          ///< The radius of the cylinder. The default value is one.

  public:

    T const & half_height() const;
    T       & half_height();
    T const & radius() const;
    T       & radius();

  public:

    Cylinder();

  public:
    EigenVector3<T> get_support_point(EigenVector3<T> dir) const override;

    T get_scale() const;

  };

  /**
   * The capsule axis is by default equal to the z-axis.
   */
  template<typename T>
  class Capsule
  : public geometry::SupportMapping<T>
  {
  public:;

  protected:

    T m_half_height;     ///< The half height of the Capsule. Default value is one.
    T m_radius;          ///< The radius of the Capsule. The default value is one.

  public:

    T const & half_height() const;
    T       & half_height();
    T const & radius()      const;
    T       & radius();

  public:

    Capsule();

  public:
      EigenVector3<T> get_support_point(EigenVector3<T> dir) const override;

    T get_scale() const;

  };

  template<typename T>
  class Ellipsoid
  : public geometry::SupportMapping<T>
  {
  protected:

    EigenVector3<T> m_scale; /// The scaling of the unit sphere along the x, y and z axes.

  public:

    void setScale(const EigenVector3<T>& vec);
    const EigenVector3<T>& scale() const;
    const EigenVector3<T>& scale();
    EigenVector3<T>& scaleRef();

public:

    Ellipsoid();

  public:
    EigenVector3<T> get_support_point(EigenVector3<T> dir) const override;

    T get_scale() const;

  };

  /**
   * A Cone.
   * A Cone is placed with its apex on positive z-axis and
   * the base orthogonal to the negative z-axis.
   *
   * The apex and base are placed equidistant wrt. x-y plane and
   * the distance is given by the half-height of the cone.
   *
   * The base is a circular disk of a specified radius. From the radius
   * and height of the cone one can compute the cone angle at the apex.
   * This angle is denoted alpha.
   */
  template<typename T>
  class Cone
  : public geometry::SupportMapping<T>
  {
  protected:

    T    m_half_height;   ///< The half height of the cone. Default value is one.
    T    m_base_radius;   ///< The radius of the circular base of the cone. Default value is one.

  public:

    T const & half_height() const;
    T       & half_height();
    T const & base_radius() const;
    T       & base_radius();

  public:

    Cone();

  public:
      EigenVector3<T> get_support_point(EigenVector3<T> dir) const override;

    T get_scale() const;

  };

  template<typename T>
  class ConvexHull
  : public geometry::SupportMapping<T>
  {

  protected:

      std::vector<EigenVector3<T>> m_points;

  public:
    void add_point(const EigenVector3<T>& p);


    const EigenVector3<T>& get_point( std::size_t const & idx) const;

    std::size_t size() const;

    void clear();

  public:

    ConvexHull();

  public:
      EigenVector3<T> get_support_point(EigenVector3<T> dir) const override;

    T get_scale() const;

  };

} // namespace convex

// CONVEX_SHAPES_H
#endif
