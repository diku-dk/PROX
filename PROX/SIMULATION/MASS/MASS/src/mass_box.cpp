#include <mass.h>

#include <tiny_is_finite.h>
#include <tiny_is_number.h>
#include <tiny_value_traits.h>

#include <cassert>

namespace mass
{

  template<typename T>
  Properties<T> compute_box(
                   T const & density
                   , T const & half_width
                   , T const & half_height
                   , T const & half_depth
                   )
  {
    using namespace tiny;

    assert( is_number( density )                || !"density must be a number"           );
    assert( is_finite( density )                || !"density must be a finite number"    );
    assert( density > 0      || !"density must be positive"           );
    assert( is_number( half_width )             || !"half_width must be a  number"       );
    assert( is_finite( half_width )             || !"half_width must be a finite number" );
    assert( half_width >= 0  || !"half_width must be non-negative"    );
    assert( is_number( half_height )            || !"half_height must be a number"       );
    assert( is_finite( half_height )            || !"half_height must be a finite number");
    assert( half_height >= 0 || !"half_height must be non-negative"   );
    assert( is_number( half_depth )             || !"half_depth must be a number"        );
    assert( is_finite( half_depth )             || !"half_depth must be a finite number" );
    assert( half_depth >= 0  || !"half_depth must be non-negative"    );

    const auto h0   = half_width;
    const auto h1   = half_height;
    const auto h2   = half_depth;
    const auto volume = 8*h0*h1*h2;
    const auto mass   = density*volume;

    Properties<T> value;

    value.m_m    = mass;
    value.m_Ixx  = ( (h1*h1 + h2*h2 ) * mass) / 3;
    value.m_Iyy  = ( (h0*h0 + h2*h2 ) * mass) / 3;
    value.m_Izz  = ( (h0*h0 + h1*h1 ) * mass) / 3;

    return value;
  }

  template Properties<float> compute_box<float>(
                                    float const & density
                                    , float const & half_width
                                    , float const & half_height
                                    , float const & half_depth
                                    );

  template Properties<double> compute_box<double>(
                                     double const & density
                                     , double const & half_width
                                     , double const & half_height
                                     , double const & half_depth
                                     );

} // namespace mass
