#include <big_stagnation.h>

#include <big_is_number.h>

#include <cmath>
#include <cassert>

namespace big
{
  
  template < typename T >
  bool stagnation(
                  ublas::vector<T> const & x_old
                  , ublas::vector<T> const & x
                  , T const & tolerance
                  )
  {
    using std::fabs;
    using std::max;
    
    
    assert( tolerance >= 0.0         || !"stagnation(): tolerance must be non-negative");
    assert( is_number( tolerance )   || !"stagnation(): internal error, NAN is encountered?");
    
    size_t const m = x.size();
    
    assert( m>0                      || !"stagnation(): zero size of x");
    assert( x.size() == x_old.size() || !"stagnation(): incompatible dimensions of x and x_old");
    
    T max_dx = 0.0;
    
    for (size_t i = 0u; i < m; ++ i)
      max_dx = max( max_dx, fabs( x(i) - x_old(i) ) );
    
    assert( is_number( max_dx ) || !"stagnation(): internal error, NAN is encountered?");
    
    if(max_dx <= tolerance)
      return true;
    
    return false;
  }
  
  
  template
  bool stagnation<float>(
                         ublas::vector<float> const & x_old
                         , ublas::vector<float> const & x
                         , float const & tolerance
                         );
  
  template
  bool stagnation<double>(
                          ublas::vector<double> const & x_old
                          , ublas::vector<double> const & x
                          , double const & tolerance
                          );
  
  
} // namespace big
