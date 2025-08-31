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

template<typename Number>
inline EigenVector3<Number> rotate(const EigenQuaternion<Number>& Quat, const EigenVector3<Number>& rVec)
{
    //return prod(  prod(q , r)  , conj(q) ).imag();
    return Quat*rVec;
}

template<typename Number>
static EigenQuaternion<Number> Rotateu (const Number rad, const EigenVector3<Number>& axis )
{
    using std::cos;
    using std::sin;
    const Number theta = rad/2;
    const Number ctheta = ( cos(theta) );
    const Number stheta = ( sin(theta) );
    EigenVector3<Number> vec = (axis).normalized() * stheta;
    return EigenQuaternion<Number>( ctheta, vec.x(), vec.y(), vec.z());
}

template<typename Number>
static EigenQuaternion<Number> Rotatex ( const Number rad)
{
    using std::cos;
    using std::sin;
    const Number theta = rad/2;
    const Number ctheta = ( cos(theta) );
    const Number stheta = ( sin(theta) );
    return EigenQuaternion<Number>(ctheta, stheta, 0, 0);
}

template<typename Number>
static EigenQuaternion<Number> Rotatey ( const Number rad)
{
    using std::cos;
    using std::sin;
    const Number theta = rad/2;
    const Number ctheta =( cos(theta) );
    const Number stheta = ( sin(theta) );
    return EigenQuaternion<Number>(ctheta, 0, stheta, 0);
}

template<typename Number>
static EigenQuaternion<Number> Rotatez ( const Number rad)
{
    using std::cos;
    using std::sin;
    const Number theta = rad/2;
    const Number ctheta = numeric_cast( cos(theta) );
    const Number stheta = numeric_cast( sin(theta) );
    return EigenQuaternion<Number>(ctheta, 0, 0, stheta);
}

/*template<typename Number>
inline EigenVector3<Number> rotate(const EigenQuaternion<Number>& Quat, const Eigen::Matrix<Number, 3,1>::NeagtiveReturnType& rVec)
{
    //return prod(  prod(q , r)  , conj(q) ).imag();
    return Quat*rVec;
}*/

template <typename TypeParameter>
inline auto coordSysToEigen(const tiny::CoordSys<TypeParameter>& input)
{
    using Number = std::remove_cv_t<std::remove_reference_t<decltype(input.T()[0])>>;
    return CoordSysEigen<Number>{toEigen(input.T()), toEigen(input.Q())};
}


template<typename T>
inline void getAxisAngle(const EigenQuaternion<T>& Q,EigenVector3<T>& axis, T& theta)
{
    using std::atan2;

    //
    // By definition a unit quaternion Q can be written as
    //
    //    Q = [s,v] = [cos(theta/2), n sin(theta/2)]
    //
    // where n is a unit vector. This is the same as a rotation of
    // theta radian around the axis n.
    //
    //
    // Rotations are difficult to work with for several reasons.
    //
    // Firstly both Q and -Q represent the same rotation. This is
    // easily proven, rotate a arbitary vector r by Q then we have
    //
    //   r^\prime = Q r Q^*
    //
    // Now rotate the same vector by -Q
    //
    //   r^\prime = (-Q) r (-Q)^* = Q r Q^*
    //
    // because -Q = [-s,-v] and (-Q)^* = [-s , v] = - [s,-v]^* = - Q^*.
    //
    // Thus the quaternion representation of a single rotation is not unique.
    //
    // Secondly the rotation it self is not well-posed. A rotation of theta
    // radians around the unit axis n could equally well be done as a rotation
    // of -theta radians around the negative unit axis n.
    //
    // This is seen by straightforward substitution
    //
    //  [ cos(-theta/2), sin(-theta/2) (-n) ] = [ cos(theta/2), sin(theta/2) n ]
    //
    // Thus we get the same quaternion regardless of whether we
    // use (+theta,+n) or (-theta,-n).
    //
    //
    // From the Quaternion we see that
    //
    //   \frac{v}{\norm{v}}  = \frac{ sin(theta/2) n }{| sin(theta/2) | } = sign(sin(theta/2)) n
    //
    // Thus we can easily get the rotation axis. However, we can not immediately
    // determine the positive rotation axis direction. The problem boils down to the
    // fact that we can not see the sign of the sinus-factor.
    //
    // Let us proceed by setting
    //
    //   x =    cos(theta/2)   =  s
    //   y =  | sin(theta/2) | =  \norm{v}
    //
    // Then we basically have two possibilities for finding theta
    //
    //  theta_1 = 2 atan2( y, x)        equivalent to      sign(sin(theta/2)) = 1
    //
    // or
    //
    //  theta_2 = 2 atan2( -y, x)       equivalent to      sign(sin(theta/2)) = -1
    //
    // If theta_1 is the solution we have
    //
    //  n = \frac{v}{\norm{v}}
    //
    // If theta_2 is the solution we must have
    //
    //  n = - \frac{v}{\norm{v}}
    //
    // Observe that we always have theta_2 = 2 pi - theta_1. Therefore theta_1 < theta_2.
    //
    // Let us imagine that we always choose $theta_1$ as the solution then
    // the correspoding quaternion for that solution would be
    //
    //         Q_1 = [cos(theta_1/2),  sin(theta_1/2)   \frac{v}{\norm{v}}]
    //             = [s ,  \norm{v}   \frac{v}{\norm{v}}]
    //             = Q
    //
    // Now if we choose theta_2 as the solution we would have
    //
    //         Q_2 = [cos(theta_2/2),  sin(theta_2/2)   -\frac{v}{\norm{v}}]
    //             = [s ,  -\norm{v}   -\frac{v}{\norm{v}}]
    //             = [s ,  \norm{v}   \frac{v}{\norm{v}}]
    //             = Q
    //
    // Thus we observe that regardless of which solution we pick we always have Q = Q_1 = Q_2.
    //
    // At this point one may be confused. However, it should be clear that theta_2 is equivalent
    // to the theta_1 rotation. The difference is simply that theta_2 corresponds to flipping the
    // rotation axis of the theta_1 case.
    //
    const T ct2   = Q.w();           //---   cos(theta/2)
    const T st2   = ( EigenVector3<T>(Q.x(), Q.y(), Q.z()) ).norm();   //---  |sin(theta/2)|

    theta = 2* atan2(st2,ct2);

    assert( st2 >= 0   || !"get_axis_angle(): |sin(theta/2)| must be non-negative");
    assert( theta >= 0 || !"get_axis_angle(): theta must be non-negative");
    assert( is_number(theta)              || !"get_axis_angle(): NaN encountered");

    axis = st2 > 0 ? EigenVector3<T>(Q.x(), Q.y(), Q.z()) / st2 : EigenVector3<T>( 0,0,0 );
}

//TINY_MATH_TYPES_H
#endif
