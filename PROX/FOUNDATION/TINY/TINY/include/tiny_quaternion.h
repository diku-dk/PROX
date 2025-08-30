#ifndef TINY_QUATERNION_H
#define TINY_QUATERNION_H

#include <tiny_value_traits.h>
#include <tiny_container.h>
#include <tiny_vector.h>
#include <tiny_random.h>

#include <cmath>          // for sin and cos

namespace tiny
{
  template < typename T >
  class Quaternion
  : public detail::Container< 1, 4, T >
  {
  protected:
      using base_class_type = detail::Container<1, 4, T>;

      enum
      {
          J_padded = base_class_type::J_padded
      };

  public:
      using real_type = typename T::real_type;
      using op_type = typename T::op_type;
      using value_traits = ValueTraits<real_type>;
      using vector3_type = Vector<3, T>;

  public:
      Quaternion()
          : base_class_type()
      {}

    ~Quaternion() {}

    Quaternion( Quaternion const & q )
		: base_class_type(q)
    {}

    Quaternion & operator = (Quaternion const & q)
		{
      base_class_type::operator=(q);
		  return *this;
		}

		Quaternion(real_type const & real, real_type const & x, real_type const & y, real_type const & z)
    : base_class_type()
		{
		  this->m_data[0] = x;
		  this->m_data[1] = y;
		  this->m_data[2] = z;
		  this->m_data[3] = real;
		}

    Quaternion(real_type const & s, vector3_type const & V )
    : base_class_type()
		{
      this->m_data[0] = V(0);
		  this->m_data[1] = V(1);
		  this->m_data[2] = V(2);
		  this->m_data[3] = s;
		}

    real_type       & real (  )       { return this->m_data[3]; }
    real_type const & real (  ) const { return this->m_data[3]; }

    vector3_type & imag (  )
		{
		  return *((vector3_type*)&(this->m_data[0]));
		}

    vector3_type  imag (  ) const
    {	
      vector3_type v;

      for (size_t i=0;i<3;++i)
        v[i] = this->m_data[i];

      return v;
    }

  public:

    static Quaternion Ru ( real_type const & rad, Vector<3,T> const & axis )
    {
      using std::cos;
      using std::sin;
      real_type const theta = rad/2;
      real_type const ctheta = value_traits::numeric_cast( cos(theta) );
      real_type const stheta = value_traits::numeric_cast( sin(theta) );
      return Quaternion ( ctheta, unit(axis) * stheta );
    }

    static Quaternion Rx ( real_type const & rad)
    {
      using std::cos;
      using std::sin;
      real_type const theta = rad/2;
      real_type const ctheta = value_traits::numeric_cast( cos(theta) );
      real_type const stheta = value_traits::numeric_cast( sin(theta) );
      return Quaternion(ctheta, stheta, 0, 0);
    }

    static Quaternion Ry ( real_type const & rad)
    {
      using std::cos;
      using std::sin;
      real_type const theta = rad/2;
      real_type const ctheta = value_traits::numeric_cast( cos(theta) );
      real_type const stheta = value_traits::numeric_cast( sin(theta) );
      return Quaternion(ctheta, 0, stheta, 0);
    }

    static Quaternion Rz ( real_type const & rad)
    {
      using std::cos;
      using std::sin;
      real_type const theta = rad/2;
      real_type const ctheta = value_traits::numeric_cast( cos(theta) );
      real_type const stheta = value_traits::numeric_cast( sin(theta) );
      return Quaternion(ctheta, 0, 0, stheta);
    }

    static Quaternion identity ( )
    {
      return  Quaternion(1, 0, 0, 0);
    }

    static Quaternion random ( real_type const & lower, real_type const & upper)
    {
      typedef typename Quaternion::iterator   iterator;

      Random<real_type> value(lower,upper);

      Quaternion q;

      iterator i     = q.begin();
      iterator end   = q.end();
      for (; i!=end ;++i)
        *i = value();

      return q;
    }

    static Quaternion random () { return random(0, 1);  }

  };

} // namespace tiny

// TINY_QUATERNION_H
#endif

