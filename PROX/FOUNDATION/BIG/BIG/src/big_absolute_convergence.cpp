#include <big_absolute_convergence.h>

#include <big_is_number.h>

#include <cassert>

namespace big
{
  template < typename T >
  bool absolute_convergence(
                            T const & f
                            , T const & tolerance
                            )
  {
    assert( tolerance >= 0              || !"absolute_convergence(): tolerance must be non-negative");
    assert( is_number( tolerance )      || !"absolute_convergence(): internal error, NAN is encountered?");
    assert( is_number( f )              || !"absolute_convergence(): internal error, NAN is encountered?");
    
    if(f <= tolerance)
      return true;
    
    return false;
  }
  
  template
  bool absolute_convergence<float>(
                                   float const & f
                                   , float const & tolerance
                                   );
  
  template
  bool absolute_convergence<double>(
                                    double const & f
                                    , double const & tolerance
                                    );
  
  
} // namespace big
