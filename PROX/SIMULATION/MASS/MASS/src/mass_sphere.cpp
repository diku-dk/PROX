#include <mass.h>

#include <tiny_is_finite.h>
#include <tiny_is_number.h>
#include <tiny_value_traits.h>

#include <cassert>

namespace mass
{

  template<typename T>
  Properties<T> compute_sphere( T const & density, T const & radius )
  {
    using namespace tiny;

    typedef ValueTraits<T> value_traits;

    assert( is_number( density )            || !"density must be a number"       );
    assert( is_finite( density )            || !"density must be a finite number");
    assert( density > 0  || !"density must be positive"       );
    assert( is_number( radius )             || !"radius must be a  number"       );
    assert( is_finite( radius )             || !"radius must be a finite number" );
    assert( radius > 0   || !"radius must be positive"        );

    T const volume  = (4*value_traits::pi()*radius*radius*radius)/3;
    T const mass    = density * volume;
    T const five    = value_traits::numeric_cast( 5.0 );

    Properties<T> value;

    value.m_m    = mass;
    value.m_Ixx  = (2*mass*radius*radius)/five;
    value.m_Iyy  = value.m_Ixx;
    value.m_Izz  = value.m_Ixx;

    return value;
  }

  template Properties<float> compute_sphere<float>(float const & density, float const & radius);

  template Properties<double> compute_sphere<double>( double const & density, double const & radius );

} // namespace mass
