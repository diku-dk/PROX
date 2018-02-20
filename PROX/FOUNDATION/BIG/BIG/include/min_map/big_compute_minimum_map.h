#ifndef BIG_COMPUTE_MINIMUM_MAP_H
#define BIG_COMPUTE_MINIMUM_MAP_H

#include <big_types.h>

namespace big
{
  
  /**
   * Minimum Map Reformulation.
   *
   * @param y    This vector contains the value of the function y = f(x).
   * @param x    The current value of the x-vector.
   * @param H    Upon return this vector holds the value of H_i(x) =  min(x_i, y_i).
   */
  template < typename T >
  void compute_minimum_map(
                           ublas::vector<T> const & y
                           , ublas::vector<T> const & x
                           , ublas::vector<T> & H
                           );

} // namespace big

// BIG_COMPUTE_MINIMUM_MAP_H
#endif
