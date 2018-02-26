#ifndef BIG_COMPUTE_PARTITIONED_JACOBIAN_H
#define BIG_COMPUTE_PARTITIONED_JACOBIAN_H

#include <big_types.h>

namespace big
{
  
  /**
   * Compute Partitioned Jacobian.
   *
   * @param A             The coefficient matrix of the linear equation y = A x + b
   * @param bitmask       A bitmask.
   * @param old2new       This vector holds an index permutation
   * @param cnt_active    This argument holds the total number of variables in the set of active constraints.
   * @param cnt_inactive  This argument holds the total number of variables in the union of the set of lower and the set of upper constraints.
   * @param A_aa          Upon return this argument holds the sub-block of reordered jacobian corresponding to the row and column sets (active,active).
   * @param A_ab          Upon return this argument holds the sub-block of reordered jacobian corresponding to the row and column sets (active,inactive).
   */
  template<typename T>
  void  compute_partitioned_jacobian(
                                     ublas::compressed_matrix<T> const & A
                                     , ublas::vector<size_t> const & bitmask
                                     , ublas::vector<size_t> const & old2new
                                     , size_t const & cnt_active
                                     , size_t const & cnt_inactive
                                     , ublas::compressed_matrix<T> & A_aa
                                     , ublas::compressed_matrix<T> & A_ab
                                     );

} // namespace big

// BIG_COMPUTE_PARTITIONED_JACOBIAN_H
#endif
