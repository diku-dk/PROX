#ifndef BIG_COMPUTE_INDEX_SETS_H
#define BIG_COMPUTE_INDEX_SETS_H

#include <big_types.h>

namespace big
{
  
  /**
   * Compute index sets.
   *
   * @param y             This vector contains the value of the function y = f(x).
   * @param x             The current value of the x-vector.
   * @param bitmask       A Bitmask with active flags
   * @param cnt_active    Upon return this argument holds the total number of active variables.
   * @param cnt_inactive  Upon return this argument holds the total number of inactive variables.
   */
  template < typename T>
  void compute_index_sets(
                          ublas::vector<T> const & y
                          , ublas::vector<T> const & x
                          , ublas::vector<size_t> & bitmask
                          , size_t  & cnt_active
                          , size_t  & cnt_inactive
                          );
  
} // namespace big

// BIG_COMPUTE_INDEX_SETS_H
#endif
