#ifndef BIG_STAGNATION_H
#define BIG_STAGNATION_H

#include <big_types.h>

namespace big
{
  
  /**
   * Stagnation Test Function.
   *
   * @param x_old                 The previous iterate value.
   * @param x                     The next iterate value.
   * @param tolerance             This argument holds the value used in the stagnation test. It is
   *                              an upper bound of the infinity-norm of the difference in the x-solution
   *                              between two iterations.  Setting the value to zero will make the test in-effective.
   * @return                      If stagnation test passes then the return value is true otherwise it is false.
   */
  template < typename T >
  bool stagnation(
                  ublas::vector<T> const & x_old
                  , ublas::vector<T> const & x
                  , T const & tolerance
                  );

} // namespace big

// BIG_STAGNATION_H
#endif
