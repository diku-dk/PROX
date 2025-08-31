#ifndef TINY_MATH_TYPES_H
#define TINY_MATH_TYPES_H

#include <tiny_vector.h>
#include <tiny_matrix.h>
#include <tiny_quaternion.h>
#include <tiny_coordsys.h>
#include <tiny_value_traits.h>
#include <tiny_type_traits.h>

#include <type_traits>

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

template <typename Number>
using EigenQuaternion = Eigen::Quaternion<Number>;


template <typename Number>
using EigenMatrix3 = Eigen::Matrix<Number, 3, 3>;

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

inline EigenQuaternion<float> toEigen(typename tiny::MathTypes<float>::quaternion_type input)
{
    typename tiny::MathTypes<float>::vector3_type vec;
    return {input.real(), vec(0), vec(1), vec(2)};
}

inline EigenQuaternion<double> toEigen(typename tiny::MathTypes<double>::quaternion_type input)
{
    typename tiny::MathTypes<double>::vector3_type vec;
    return {input.real(), vec(0), vec(1), vec(2)};
}

inline typename tiny::MathTypes<float>::quaternion_type fromEigen(EigenQuaternion<float> input)
{
    return tiny::MathTypes<float>::quaternion_type(input.w(), input.x(), input.y(), input.z());
}

inline typename tiny::MathTypes<double>::quaternion_type fromEigen(EigenQuaternion<double> input)
{
    return tiny::MathTypes<double>::quaternion_type(input.w(), input.x(), input.y(), input.z());
}



template <typename Vector>
inline auto dot(Vector a, Vector b)
{
    return a.dot(b);
}

template< typename Number >
class CoordSysEigen
{
public:

protected:
    EigenVector3<Number> m_T;      ///< The Position.
    EigenQuaternion<Number> m_Q;      ///< The orientation in Quaternion form.

public:

    EigenVector3<Number>& T()       { return m_T; }
    const EigenVector3<Number>& T() const { return m_T; }

    EigenQuaternion<Number>& Q()       { return m_Q; }
    const EigenQuaternion<Number>& Q() const { return m_Q; }

public:

    CoordSysEigen()
        : m_T( 0 )
        , m_Q( 1, 0, 0, 0)
    {}

    CoordSysEigen(const CoordSysEigen& X)
        : m_T(X.m_T)
        , m_Q(X.m_Q)
    { }


    CoordSysEigen(EigenVector3<Number> const & T_val,const EigenQuaternion<Number>& Q_val)
    {
        m_T = T_val;
        m_Q = (Q_val).normalized();
    }

    CoordSysEigen(EigenVector3<Number> const & T_val, const EigenMatrix3<Number>& R_val)
    {
        m_T = T_val;
        EigenQuaternion<Number> Quat(R_val);
        m_Q = Quat;
    }

    CoordSysEigen & operator=(CoordSysEigen const & C)
    {
        m_T = C.m_T;
        m_Q = C.m_Q;
        return *this;
    }

    bool operator==(CoordSysEigen const & C) const {  return m_T == C.m_T && m_Q==C.m_Q; }

public:

    static CoordSysEigen identity()
    {
        CoordSysEigen X;
        X.T() = {0,0,0};
        X.Q() = EigenQuaternion<Number>::Identity();
        return X;
    }

    static CoordSysEigen make(const EigenVector3<Number> & T, const EigenQuaternion<Number>& Q)
    {
        CoordSysEigen X;
        X.T() = T;
        X.Q() = Q;
        return X;
    }

};

template <typename TypeParameter>
inline auto coordSysToEigen(const tiny::CoordSys<TypeParameter>& input)
{
    using Number = std::remove_cv_t<std::remove_reference_t<decltype(input.T()[0])>>;
    return CoordSysEigen<Number>{toEigen(input.T()), toEigen(input.Q())};
}


//TINY_MATH_TYPES_H
#endif
