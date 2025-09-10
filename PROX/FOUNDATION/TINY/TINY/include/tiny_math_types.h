#ifndef TINY_MATH_TYPES_H
#define TINY_MATH_TYPES_H

#include <eigen3/Eigen/Dense>
#include <numbers>
#include <random>
#include <tiny_coordsys.h>
#include <tiny_matrix.h>
#include <tiny_quaternion.h>
#include <tiny_type_traits.h>
#include <tiny_value_traits.h>
#include <tiny_vector.h>
#include <type_traits>

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
#endif // TINY_MATH_TYPES_H
