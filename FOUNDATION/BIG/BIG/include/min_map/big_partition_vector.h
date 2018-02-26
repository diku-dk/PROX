#ifndef BIG_PARTITION_VECTOR_H
#define BIG_PARTITION_VECTOR_H

#include <big_types.h>

namespace big
{
  /**
   * Partition vector.
   *
   * @param x             The value of the vector.
   * @param bitmask       A bitmask of active flags.
   * @param old2new       This vector holds an index permutation of old indices
   *                      to new indices. The permutation is such that the first sub-block of
   *                      the new vector corresponds to active constraints. The second corresponds
   *                      to lower constraints and the last sub-block corresponds to the upper
   *                      constraints set.
   * @param cnt_active    This argument holds the total number of active variables.
   * @param cnt_inactive  This argument holds the total number of inactive variables.
   * @param x_a           Upon return this vector holds the values of the sub-block of rhs that corresponds to the active constraints.
   * @param x_b           Upon return this vector holds the values of the sub-block of rhs that corresponds to the inactive constraints.
   */
  template<typename T >
  void  partition_vector(
                         ublas::vector<T> const & x
                         , ublas::vector<size_t> const & bitmask
                         , ublas::vector<size_t> const & old2new
                         , size_t const & cnt_active
                         , size_t const & cnt_inactive
                         , ublas::vector<T>  & x_a
                         , ublas::vector<T>  & x_b
                         );

} // namespace big

// BIG_PARTITION_VECTOR_H
#endif
