#include <big_relative_convergence.h>

#include <big_is_number.h>

#include <cmath>
#include <cassert>

namespace big
{
  
  template < typename T >
  bool relative_convergence(
                            T const & f_old
                            , T const & f
                            , T const & tolerance
                            )
  {
    using std::fabs;
    
    assert( tolerance >= 0            || !"relative_convergence(): tolerance must be non-negative");
    assert( is_number( tolerance )    || !"relative_convergence(): internal error, NAN is encountered?");
    assert( is_number( f_old )        || !"relative_convergence(): internal error, NAN is encountered?");
    assert( is_number( f )            || !"relative_convergence(): internal error, NAN is encountered?");
    
    
    if(f_old == 0)
    {
      // f and f_old are both zero
      if(f == 0)
        return true;
      // f_old is zero and f is non-zero, so the test does not really make sense!
      return false;
    }
    
    T const relative_test = fabs(f - f_old) / fabs(f_old);
    
    assert( is_number( relative_test ) || !"relative_convergence(): internal error, NAN is encountered?");
    
    if(relative_test <= tolerance)
      return true;
    
    return false;
  }
  
  template
  bool relative_convergence<float>(
                                   float const & f_old
                                   , float const & f
                                   , float const & tolerance
                                   );
  
  template
  bool relative_convergence<double>(
                                    double const & f_old
                                    , double const & f
                                    , double const & tolerance
                                    );
  
  
  
} // namespace big
