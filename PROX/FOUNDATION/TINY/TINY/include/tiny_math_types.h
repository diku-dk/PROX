#ifndef TINY_MATH_TYPES_H
#define TINY_MATH_TYPES_H

#include <tiny_vector.h>
#include <tiny_matrix.h>
#include <tiny_quaternion.h>
#include <tiny_coordsys.h>
#include <tiny_value_traits.h>
#include <tiny_type_traits.h>

#include <eigen3/Eigen/Dense>

namespace tiny
{

  /**
   * Math Types Container.
   * This class contains the most common used small-sized matrix
   * types. It is intended as a convience tool making. It makes it easier to
   * pass math types around as a small collection of types for instance as
   * template arguments for data structures or algorithms.
   *
   * @tparam T   The precision of the math types. Typical either float or double.
   */
  template< typename T >
  class MathTypes;

  template<  >
  class MathTypes<float>
  {
  public:
      using type_traits = float_traits;
      using real_type = type_traits::real_type;
      using vector3_type = Vector<3, type_traits>;
      using vector4_type = Eigen::Matrix<float, 4, 1>;

      using quaternion_type = Quaternion<type_traits>;
      using matrix3x3_type = Matrix<3, 3, type_traits>;
      using matrix4x4_type = Matrix<4, 4, type_traits>;
      using coordsys_type = CoordSys<type_traits>;
      using value_traits = ValueTraits<real_type>;
  };

  template<  >
  class MathTypes<double>
  {
  public:
      using type_traits = double_traits;
      using real_type = type_traits::real_type;
      using vector3_type = Vector<3, type_traits>;
      using vector4_type = Eigen::Matrix<double, 4, 1>;

      using quaternion_type = Quaternion<type_traits>;
      using matrix3x3_type = Matrix<3, 3, type_traits>;
      using coordsys_type = CoordSys<type_traits>;
      using value_traits = ValueTraits<real_type>;
  };

} // namespace tiny

template <typename Number>
using EigenVector3 = Eigen::Vector<Number, 3>;

inline EigenVector3<float> toEigen(typename tiny::MathTypes<float>::vector3_type input)
{
    return {input(0), input(1), input(2)};
}

inline EigenVector3<double> toEigen(typename tiny::MathTypes<double>::vector3_type input)
{
    return {input(0), input(1), input(2)};
}

inline typename tiny::MathTypes<float>::vector3_type fromEigen(EigenVector3<float> input)
{
    return tiny::MathTypes<float>::vector3_type::make(input.x(), input.y(), input.z());
}

inline typename tiny::MathTypes<double>::vector3_type fromEigen(EigenVector3<double> input)
{
    return tiny::MathTypes<double>::vector3_type::make(input.x(), input.y(), input.z());
}

template <typename Vector>
inline auto dot(Vector a, Vector b)
{
    return a.dot(b);
}

//TINY_MATH_TYPES_H
#endif
