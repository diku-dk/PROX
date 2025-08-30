#ifndef TINY_COORDSYS_H
#define TINY_COORDSYS_H

#include <tiny_vector.h>
#include <tiny_matrix.h>
#include <tiny_quaternion.h>
#include <tiny_quaternion_functions.h>  // for unit and make
#include <tiny_value_traits.h>

namespace tiny
{

  /**
   * A Coordinate System.
   * This class is used to represents a transform that brings you from a local frame
   * into a global frame. We write this mathematically as:
   *
   *  X : BF -> WCS
   *
   */
  template< typename T_ >
  class CoordSys
    {
    public:
        using real_type = typename T_::real_type;
        using value_traits = ValueTraits<real_type>;
        using vector3_type = Vector<3, T_>;
        using quaternion_type = Quaternion<T_>;
        using matrix3x3_type = Matrix<3, 3, T_>;

    protected:
        vector3_type m_T;      ///< The Position.
        quaternion_type m_Q;      ///< The orientation in Quaternion form.

    public:

      vector3_type       & T()       { return m_T; }
      vector3_type const & T() const { return m_T; }

      quaternion_type       & Q()       { return m_Q; }
      quaternion_type const & Q() const { return m_Q; }

    public:

      CoordSys()
      : m_T( 0 )
      , m_Q( 1, 0, 0, 0)
      {}

      CoordSys( CoordSys const & X)
      : m_T(X.m_T)
      , m_Q(X.m_Q)
      { }


      CoordSys(vector3_type const & T_val, quaternion_type const & Q_val)
      {
        m_T = T_val;
        m_Q = unit(Q_val);
      }

      CoordSys(vector3_type const & T_val, matrix3x3_type const & R_val)
      {
        m_T = T_val;
        m_Q = tiny::make( R_val );
      }

      CoordSys & operator=(CoordSys const & C)
      {
        m_T = C.m_T;
        m_Q = C.m_Q;
        return *this;
      }

      bool operator==(CoordSys const & C) const {  return m_T == C.m_T && m_Q==C.m_Q; }

    public:

      static CoordSys identity()
      {
        CoordSys X;
        X.T().zero();
        X.Q() = quaternion_type::identity();
        return X;
      }

      static CoordSys make(vector3_type const & T, quaternion_type const & Q)
      {
        CoordSys X;
        X.T() = T;
        X.Q() = Q;
        return X;
      }

    };

} // namespace tiny

//TINY_COORDSYS_H
#endif
