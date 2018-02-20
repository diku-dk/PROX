#ifndef BIG_AGGLOMERATE_VECTOR_H
#define BIG_AGGLOMERATE_VECTOR_H

#include <big_types.h>

namespace big
{
  
  /**
   * Transfer partitioned vector back to the global un-partitioned system.
   *
   * @param x_a          A sub-vector that should be transfered into the higher-dimensional vector x.
   * @param x_b          A sub-vector that should be transfered into the higher-dimensional vector x.
   *
   * @param new2old       This vector holds an index permutation of the new indices (in the sub-vectors)
   *                      to the old indices (in the un-partitioned system).
   *
   * @param x             Upon return this argument holds the agglomerated vector x = pi( x_a, x_b)
   */
  template<typename T>
  void agglomerate_vector(
                          ublas::vector<T> const & x_a
                          , ublas::vector<T> const & x_b
                          , ublas::vector<size_t> const & new2old
                          , ublas::vector<T> & x
                          );
  
} // namespace big

// BIG_AGGLOMERATE_VECTOR_H
#endif
